%% === 最终修正版：一维水平扫描 (适配你的函数定义) ===
clc; clear; close all;

% [关键] 自动添加路径
addpath(genpath(pwd)); 

% ================= 1. 系统参数 =================
ht_BS = 40; h_RIS = 3; hr_User = 8;
bs_pos = [0, 0, ht_BS];
user_X_center = 10000; user_Y_center = 1000; user_R = 200;

% 通信参数
f_c = 1.9; lambda = 3e8/(f_c*1e9); K_rician = 10;
K = 4; M = 8; N = 256; 
Ps_max = 20000; Pr_max = 500;
sigma2 = 1e-10; sigmar2 = 1e-10; eta_k = ones(K,1);

% ================= 2. 生成固定用户 =================
fprintf('>>> 正在生成用户 (固定种子)...\n');
rng(888); 
theta = 2*pi*rand(K,1);
rad   = user_R * rand(K,1);
user_pos = zeros(K,3);
for k = 1:K
    user_pos(k,:) = [user_X_center + rad(k)*cos(theta(k)), ...
                     user_Y_center + rad(k)*sin(theta(k)), ...
                     hr_User];
end

% ================= 3. 扫描参数 =================
Scan_X_Vec = 1000 : 200 : 11000; % 步长 200m
Scan_Y_Fixed = user_Y_center;
WSR_Curve = zeros(length(Scan_X_Vec), 1);

% ================= 4. 开始扫描 =================
fprintf('>>> 开始扫描 %d 个点...\n', length(Scan_X_Vec));

for i = 1:length(Scan_X_Vec)
    curr_ris_x = Scan_X_Vec(i);
    
    % [核心修正点] 去掉了最后一个 false 参数，只保留 12 个参数
    % 对应定义: build_slot_channels(K, N_RIS, M, bs_pos, user_pos, x_traj, y_traj, h_const, lambda, gamma_reflect, K_rician, no_direct_link)
    try
        [channels, ~, ~] = build_slot_channels(K, N, M, bs_pos, user_pos, ...
            curr_ris_x, Scan_Y_Fixed, h_RIS, lambda, 0, K_rician, false);
        ch = channels(1);
    catch ME
        error(['生成信道失败！\n错误原因: %s\n' ...
               '请确认 build_slot_channels.m 的参数个数是否匹配。'], ME.message);
    end
    
    % 初始化
    Theta0 = diag(exp(1j*2*pi*rand(N,1)));
    W0     = ones(K*M,1) * sqrt(Ps_max/K/M);
    
    % 计算 WSR
    try
        [~, ~, wsr_val] = active_RIS_precoding(M,K,N,Ps_max,Pr_max,sigma2,sigmar2,eta_k, Theta0, W0, ch.h_k, ch.f_k, ch.G);
        WSR_Curve(i) = wsr_val;
    catch ME
        warning('点 X=%.0f 计算失败: %s', curr_ris_x, ME.message);
        WSR_Curve(i) = 0;
    end
    
    if mod(i, 5) == 0, fprintf('进度: %.0f%%\n', i/length(Scan_X_Vec)*100); end
end

% ================= 5. 绘图 =================
figure('Color','w'); 
plot(Scan_X_Vec, WSR_Curve, 'r-o', 'LineWidth', 2);
xline(user_X_center, 'k--', 'User Center');
xlabel('RIS X Position (m)'); ylabel('Sum Rate'); grid on;
title('Deployment Scan (Corrected Version)');