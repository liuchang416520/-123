clear; clc; close all;

% ==========================================================
% 有源 RIS 1D 静态扫描测试 (寻找最优部署 X 坐标)
% ==========================================================

% 1. 基本场景参数 (与 main_fp_sca_trajectory.m 同步)
ht_BS   = 35;          % 基站高度 35m
h_RIS   = 15;          % USV (RIS) 高度 15m
hr_User = 10;          % 用户高度 10m

user_X_center = 45000; % 用户群中心 45km
user_Y_center = 0;
user_R        = 1500;

K = 4;                 % 用户数
M = 10;                % 基站天线数
N = 256;               % RIS 单元数

% Ps_max_dBm = 47;
% Ps_max = 10^(Ps_max_dBm/10);
% Pr_max = Ps_max * 0.8; % 有源 RIS 功率
P_total_dBm = 47;
P_total_mW  = 10^(P_total_dBm/10); % 总功率约束 50000 mW
RIS_power_ratio = 0.001; 

Pr_max = P_total_mW * RIS_power_ratio;      % 有源 RIS 最大发射功率
Ps_max = P_total_mW * (1 - RIS_power_ratio); % 基站最大发射功率

fprintf('总功率约束: %.2f W. 基站分配: %.2f W, RIS分配: %.2f W\n', P_total_mW/1000, Ps_max/1000, Pr_max/1000);

sigma2  = 1e-13;       % 噪声功率
sigmar2 = 1e-11;       
eta_k   = ones(K,1);   

f_c   = 5;             
lambda = 3e8/(f_c*1e9);

% 【修改1】：将 K_rician 改回 1e8。我们要看纯理论距离衰减趋势，10 会引入多径散射的微小随机起伏
K_rician = 10;        

tau_power   = 4;  
kappa_power = 4;  
gamma_reflect = 0;
bs_pos = [0, 0, ht_BS];

% 2. 生成固定的用户位置
rng(42); 
theta = 2*pi*rand(K,1);
rad   = user_R*rand(K,1);
user_pos = zeros(K,3);
for k = 1:K
    user_pos(k,:) = [user_X_center + rad(k)*cos(theta(k)), user_Y_center + rad(k)*sin(theta(k)), hr_User];
end

% 3. 扫描点设置
X_scan = 1000:1000:44000; 
num_points = length(X_scan);
WSR_scan = zeros(1, num_points);

h_k0 = zeros(M,K); f_k0 = zeros(N,K); G0 = zeros(N,M);

fprintf('=== 开始 1D 静态扫描测试 ===\n');
fprintf('总计 %d 个扫描点，请耐心等待...\n', num_points);

% ==========================================
% 【修改2】：初始化确实要在循环外部！且后续不重置
% ==========================================
rng(123); 
Theta0 = diag(exp(1j*2*pi*rand(N,1)));
W0 = exp(1j*2*pi*rand(K*M,1))*sqrt(Ps_max/K/M);

% ==========================================
% 分离缩放参数
% ==========================================
scale_G = 1e3;  
scale_f = 1e6;  

fprintf('=== 开始 1D 静态扫描测试 (终极版：满载放大 + 智能对齐 + 热启动) ===\n');

for i = 1:num_points
    x_current = X_scan(i);
    y_current = 0; 
    ris_pos = [x_current, y_current, h_RIS];
    
    d_BR = max(norm(bs_pos - ris_pos), 1.0);
    d_RU = zeros(K,1); d_BU = zeros(K,1);
    for k = 1:K
        d_RU(k) = max(norm(ris_pos - user_pos(k,:)), 1.0);
        d_BU(k) = max(norm(bs_pos - user_pos(k,:)), 1.0);
    end
    
    rng(999); % 冻结空间相位
    [h_k, f_k, G, ~, ~, ~] = Channel_generate2(K, N, M, d_BR, d_BU, d_RU, ht_BS, h_RIS, hr_User, lambda, gamma_reflect, K_rician);
    
    % --- 应用分离缩放 ---
    G_scaled = G * scale_G;
    f_k_scaled = f_k * scale_f;
    h_k_scaled = h_k * 0.354 * (scale_G * scale_f) ; % 削弱直射
    
    Pr_max_eff = Pr_max * (scale_G^2); 
    sigmar2_eff = sigmar2 * (scale_G^2);
    sigma2_eff = sigma2 * (scale_G^2 * scale_f^2);
    
    % =========================================================
    % 【全网最强初始化：打破有源 RIS 的零陷阱】
    % =========================================================
    direct_phase = angle(h_k_scaled(1,1));
        % 1. 先用模长为1的矩阵把相位对齐
    %theta_phase = - (angle(f_k_scaled(1,:).') + angle(G_scaled(:,1)));
    theta_phase = direct_phase - angle(f_k_scaled(1,:).') - angle(G_scaled(:,1));
        Theta_tmp = diag(exp(1j * theta_phase));
        
        % 2. 让基站满功率发射 (MRT)
        W0_mat = zeros(M, K);
        for k = 1:K
           % h_eff = f_k_scaled(k,:) * Theta_tmp * G_scaled;
           h_eff = h_k_scaled(k,:) + f_k_scaled(k,:) * Theta_tmp * G_scaled;
            W0_mat(:,k) = (h_eff' / norm(h_eff)) * sqrt(Ps_max / K);
        end
        
        % 3. 【绝对核心】：计算此时 RIS 接收到的总功率 (信号 + 热噪声)
        P_rx = 0;
        for k = 1:K
            P_rx = P_rx + norm(G_scaled * W0_mat(:,k))^2;
        end
        P_rx = P_rx + N * sigmar2_eff;
        
        % 4. 强行计算 RIS 所能提供的最大物理放大倍数 a_max (通常在数万倍)
        a_max = sqrt(Pr_max_eff / P_rx);
        
        % 5. 赋予有源 RIS 真正的灵魂：不仅相干，而且满载放大！
        Theta0 = a_max * Theta_tmp;
        W0 = W0_mat(:); 
    
    % =========================================================
    
    traj_cfg = struct();
    traj_cfg.enable = false; 
    
    % 执行 BCD 优化
    [W_opt, Theta_opt, Rsum] = active_RIS_precoding(...
        M, K, N, Ps_max, Pr_max_eff, sigma2_eff, sigmar2_eff, eta_k, ...
        Theta0, W0, h_k_scaled, f_k_scaled, G_scaled, traj_cfg);
    
    WSR_scan(i) = max(Rsum);
    
   
    
    fprintf('进度 %d/%d: X位置 = %5d m, 最优 WSR = %6.4f bps/Hz\n', i, num_points, x_current, WSR_scan(i));
end

% 5. 结果绘图
figure('Position', [150, 150, 800, 500]);
plot(X_scan, WSR_scan, '-o', 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5, 'MarkerSize', 7, 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
grid on;
set(gca, 'GridAlpha', 0.4, 'FontSize', 11);
xlabel('USV (Active RIS) X 轴坐标位置 (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('系统加权和速率 WSR (bps/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('有源 RIS 最佳部署位置 1D 扫描 (纯视距, 45km场景)', 'FontSize', 14, 'FontWeight', 'bold');

% 标记 WSR 的理论最高点
[max_wsr, max_idx] = max(WSR_scan);
best_x = X_scan(max_idx);
hold on;
plot(best_x, max_wsr, 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

% 添加文字说明标签
text(best_x, max_wsr * 1.02, sprintf('最佳部署点: %d m\n最高 WSR: %.2f', best_x, max_wsr), ...
    'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');

% 标注基站与用户所在的位置基准线
xline(0, '--k', '基站 (BS)', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
xline(user_X_center, '--k', '用户群中心', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
xlim([-2000, 48000]);