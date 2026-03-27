clear; clc; close all;
% ==========================================================
% 有源 vs 无源 RIS 固定位置 (812.7, -79.8) N 敏感性对比研究
% 采用单次物理信道切片 + 双路增量热启动策略
% ==========================================================

% 固定 RIS 位置（用户指定）
x_fixed = 812.7;
y_fixed = -79.8;

% 用户位置与系统参数
ht_BS = 5; h_RIS = 2; hr_User = 1.5;
K = 4; M = 10;
P_total_dBm = 20;
P_total_mW  = 10^(P_total_dBm/10);
RIS_power_ratio = 0.01;
N_ref = 100;

% 功率分配：为了绝对公平，无源 RIS 独占全部 20dBm 功率给基站
Pr_elem_base = (P_total_mW * RIS_power_ratio) / N_ref;
Ps_max_active = P_total_mW * (1 - RIS_power_ratio); % 有源时的基站功率
Ps_max_passive = P_total_mW;                        % 无源时的基站全功率

sigma2 = 1e-13; sigmar2 = 1e-13; eta_k = ones(K,1);
f_c = 5; lambda = 3e8/(f_c*1e9); K_rician = 10;
gamma_reflect = 0;
scale_G = 1e3; scale_f = 1e6;

% 有效噪声
sigma2_eff = sigma2 * (scale_G^2 * scale_f^2);
sigmar2_eff = sigmar2 * (scale_G^2);

bs_pos = [0, 0, ht_BS];
rng(42);
user_X_center = 1000; user_Y_center = 0; user_R = 50;
theta = 2*pi*rand(K,1); rad = user_R*rand(K,1);
user_pos = zeros(K,3);
for k = 1:K
    user_pos(k,:) = [user_X_center + rad(k)*cos(theta(k)), user_Y_center + rad(k)*sin(theta(k)), hr_User];
end

% N 扫描
N_list = 50:50:700;
nN = numel(N_list);

% [新增] 分别记录无源和有源的 WSR
Rsum_passive_list = zeros(nN, 1);
Rsum_active_list = zeros(nN, 1);

fprintf('========== 有源 vs 无源 RIS N 敏感性对比 ==========\n');
fprintf('场景: 直达链路阻塞, P_total=%ddBm\n', P_total_dBm);
fprintf('有源规则: 基站功率=%.2fmW, RIS单阵子恒定功率=%.4fmW\n', Ps_max_active, Pr_elem_base);
fprintf('无源规则: 基站独占全部功率=%.2fmW\n\n', Ps_max_passive);

% 锁定一个确定的随机种子，生成唯一且包含所有阵子信息的“超级大信道”
N_max = max(N_list); 
rng(999); 
[channels_max, ~, ~] = build_slot_channels(K, N_max, M, bs_pos, user_pos, x_fixed, y_fixed, h_RIS, lambda, gamma_reflect, K_rician, true, scale_G, scale_f);
ch_max = channels_max(1);

% ==========================================================
% [核心修改] 为无源和有源分别建立独立的“最优解记忆库”
Theta_passive_opt_prev = [];
W_passive_opt_prev = [];

Theta_active_opt_prev = [];
W_active_opt_prev = [];
% ==========================================================

for ii = 1:nN
    N = N_list(ii);
    fprintf('N = %d (%d/%d) ... ', N, ii, nN);

    % 动态功率设置：每个阵子功率恒定，总功率随 N 增加
    Pr_max_N = Pr_elem_base * N; 
    Pr_max_eff = Pr_max_N * (scale_G^2);

    % 信道切片（保证物理空间绝对一致）
    ch_current_G = ch_max.G(1:N, :);       
    ch_current_f_k = ch_max.f_k(:,1:N);   
    ch_current_h_k = ch_max.h_k;           

    % ==========================================================
    % 步骤 1：无源 RIS 计算 (带有增量热启动)
    % ==========================================================
    if ii == 1
        rng(123); 
        Theta_init_pass = diag(exp(1j*2*pi*rand(N,1)));
        W_init_pass = exp(1j*2*pi*rand(K*M,1)) * sqrt(Ps_max_passive/K/M);
    else
        N_prev = N_list(ii-1);
        Theta_init_pass = diag(exp(1j*2*pi*rand(N,1))); % 新阵子给随机相位
        Theta_init_pass(1:N_prev, 1:N_prev) = Theta_passive_opt_prev; % 继承老阵子最优相位
        W_init_pass = W_passive_opt_prev; % 继承预编码
    end

    [W_passive_opt_prev, Theta_passive_opt_prev, Rsum_passive_list(ii)] = passive_RIS_precoding(M, K, N, Ps_max_passive, sigma2_eff, eta_k, Theta_init_pass, W_init_pass, ch_current_h_k, ch_current_f_k, ch_current_G);

    % ==========================================================
    % 步骤 2：有源 RIS 计算 (带有增量热启动)
    % ==========================================================
    if ii == 1
        % 对于 N 的起点：借用无源的最优相位，并暴力打针 (100倍注水) 唤醒梯度
        Theta_init_act = 100 * Theta_passive_opt_prev; 
        W_init_act = W_passive_opt_prev; 
    else
        % 对于后续的 N：坚定继承上一次【有源】的最优解！
        % 给新增加的阵子赋予微小的非零随机初值防奇异报错
        N_prev = N_list(ii-1);
        Theta_init_act = 0.1 * diag(exp(1j*2*pi*rand(N,1))); 
        Theta_init_act(1:N_prev, 1:N_prev) = Theta_active_opt_prev; 
        W_init_act = W_active_opt_prev; 
    end

    [W_active_opt_prev, Theta_active_opt_prev, Rsum_active_list(ii)] = active_RIS_precoding(M, K, N, Ps_max_active, Pr_max_eff, sigma2_eff, sigmar2_eff, eta_k, Theta_init_act, W_init_act, ch_current_h_k, ch_current_f_k, ch_current_G);

    fprintf('无源WSR=%.2f \t 有源WSR=%.2f\n', Rsum_passive_list(ii), Rsum_active_list(ii));
end

% ==========================================================
% 绘图：双子图，绝对 WSR 对比 + 相对增益
% ==========================================================
figure('Position', [50, 50, 1000, 420]);

subplot(1, 2, 1); hold on; grid on; box on;
plot(N_list, Rsum_active_list, '-r', 'LineWidth', 2, 'Marker', 's', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
plot(N_list, Rsum_passive_list, '-b', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
xlabel('反射元件数量 N'); ylabel('WSR (bps/Hz)');
title('有源 vs 无源 WSR 对比');
legend('Active RIS', 'Passive RIS', 'Location', 'best');
yr = max(Rsum_active_list) - min(Rsum_passive_list);
if yr < 0.5, yr = 0.5; end  
ylim([min(Rsum_passive_list)-0.1*yr, max(Rsum_active_list)+0.1*yr]);

subplot(1, 2, 2); hold on; grid on; box on;
wsr_ref_active = Rsum_active_list(1);
wsr_gain_pct_active = (Rsum_active_list - wsr_ref_active) / wsr_ref_active * 100;
plot(N_list, wsr_gain_pct_active, '-r', 'LineWidth', 2, 'Marker', 's', 'MarkerSize', 6, 'MarkerFaceColor', 'r');

wsr_ref_passive = Rsum_passive_list(1);
wsr_gain_pct_passive = (Rsum_passive_list - wsr_ref_passive) / wsr_ref_passive * 100;
plot(N_list, wsr_gain_pct_passive, '-b', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');

plot([N_list(1), N_list(end)], [0, 0], '--k', 'LineWidth', 1);
xlabel('反射元件数量 N'); ylabel('相对增益 (%)');
title(sprintf('相对 N=%d 的自我提升比例', N_list(1)));
legend('Active RIS 增幅', 'Passive RIS 增幅', 'Location', 'best');

save('fixed_812_N_results.mat', 'x_fixed', 'y_fixed', 'N_list', 'Rsum_active_list', 'Rsum_passive_list');
fprintf('\n========== 研究完成，结果已保存至 fixed_812_N_results.mat ==========\n');