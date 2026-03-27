clear; clc; close all;
% ==========================================================
% 有源 RIS 反射单元数量 N 敏感性研究：N = 64:2:128
% 研究轨迹和 WSR 随 N 的变化（不做无源 RIS）
% ==========================================================

% 1. 轨迹起终点
x_start = 50;    y_start = -300;
x_end   = 900;   y_end   = -300;

% 2. 边界限制（仅 y 方向）
y_min = -300; y_max = 0;

% 3. 运动学参数
v_max = 10;      T_total = 120;   Nslot = 50;
dt = T_total / Nslot;

% 4. 用户位置
user_X_center = 1000; user_Y_center = 0; user_R = 50;

% 5. 高度与通信设置（N 在循环中变化）
ht_BS = 5; h_RIS = 2; hr_User = 1.5;
K = 4; M = 10;
P_total_dBm = 15;
P_total_mW  = 10^(P_total_dBm/10);
RIS_power_ratio = 0.1;
Pr_max = P_total_mW * RIS_power_ratio;
Ps_max = P_total_mW * (1 - RIS_power_ratio);

sigma2 = 1e-13; sigmar2 = 1e-13; eta_k = ones(K,1);
f_c = 5; lambda = 3e8/(f_c*1e9); K_rician = 10;
tau_power = 4; kappa_power = 4; gamma_reflect = 0;
scale_G = 1e3; scale_f = 1e6;

% 初始轨迹
x_traj = linspace(x_start, x_end, Nslot).';
y_traj = linspace(y_start, y_end, Nslot).';

bs_pos = [0, 0, ht_BS];
rng(42);
theta = 2*pi*rand(K,1); rad = user_R*rand(K,1);
user_pos = zeros(K,3);
for k = 1:K
    user_pos(k,:) = [user_X_center + rad(k)*cos(theta(k)), user_Y_center + rad(k)*sin(theta(k)), hr_User];
end

% N 扫描：64, 66, ..., 128
N_list = 64:16:128;
nN = numel(N_list);
Rsum_list = zeros(nN, 1);
wsr_mean_list = zeros(nN, 1);
wsr_max_list = zeros(nN, 1);
x_opt_cell = cell(nN, 1);
y_opt_cell = cell(nN, 1);
rate_hist_cell = cell(nN, 1);

fprintf('========== 有源 RIS N 敏感性研究 (N=%d:%d:%d) ==========\n', N_list(1), 2, N_list(end));

for ii = 1:nN
    N = N_list(ii);
    fprintf('\n--- N = %d (%d/%d) ---\n', N, ii, nN);

    % 信道维度随 N 变化：h_k M×K, f_k N×K, G N×M
    rng(123);
    h_k0 = zeros(M, K);
    f_k0 = zeros(N, K);
    G0 = zeros(N, M);
    Theta0 = diag(exp(1j*2*pi*rand(N,1)));
    W0 = exp(1j*2*pi*rand(K*M,1)) * sqrt(Ps_max/K/M);

    traj_cfg = struct();
    traj_cfg.enable = true;
    traj_cfg.bs_pos = bs_pos; traj_cfg.user_pos = user_pos;
    traj_cfg.x_traj = x_traj; traj_cfg.y_traj = y_traj;
    traj_cfg.h_const = h_RIS; traj_cfg.Nslot = Nslot;
    traj_cfg.dt = dt; traj_cfg.v_max = v_max;
    traj_cfg.y_min = y_min; traj_cfg.y_max = y_max;
    traj_cfg.lambda = lambda; traj_cfg.gamma_reflect = gamma_reflect;
    traj_cfg.K_rician = K_rician; traj_cfg.tau_power = tau_power; traj_cfg.kappa_power = kappa_power;
    traj_cfg.no_direct_link = false;
    traj_cfg.MaxIter = 15; traj_cfg.MinIter = 5;
    traj_cfg.return_to_endpoint = true;
    traj_cfg.scale_G = scale_G; traj_cfg.scale_f = scale_f;
    traj_cfg.no_plot = true;  % 批量研究时关闭单次绘图

    [~, ~, Rsum, ~, ~, x_opt, y_opt, rate_hist, ~] = active_RIS_precoding(M, K, N, Ps_max, Pr_max, sigma2, sigmar2, eta_k, Theta0, W0, h_k0, f_k0, G0, traj_cfg);

    wsr_final = rate_hist(:, end);
    Rsum_list(ii) = Rsum;
    wsr_mean_list(ii) = mean(wsr_final);
    wsr_max_list(ii) = max(wsr_final);
    x_opt_cell{ii} = x_opt;
    y_opt_cell{ii} = y_opt;
    rate_hist_cell{ii} = wsr_final;

    fprintf('  Rsum=%.4f, Mean WSR=%.4f, Max WSR=%.4f\n', Rsum, wsr_mean_list(ii), wsr_max_list(ii));
end

% ==========================================================
% 结果绘图
% ==========================================================
figure('Position', [50, 50, 1200, 500]);

% 子图一：WSR vs N
subplot(1, 2, 1); hold on; grid on; box on;
plot(N_list, wsr_mean_list, '-b', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
plot(N_list, wsr_max_list, '--r', 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
xlabel('反射单元数量 N');
ylabel('WSR (bps/Hz)');
legend({'平均 WSR', '最大 WSR'}, 'Location', 'best');
title('有源 RIS：WSR 随 N 变化');

% 子图二：轨迹对比（选取若干代表性 N）
subplot(1, 2, 2); hold on; grid on; box on;
plot(bs_pos(1), bs_pos(2), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
scatter(user_pos(:,1), user_pos(:,2), 50, 'b', 'filled');
th = linspace(0, 2*pi, 100);
plot(user_X_center + user_R*cos(th), user_Y_center + user_R*sin(th), '--', 'Color', [0.6 0.6 0.6]);
plot(x_traj, y_traj, '--k', 'LineWidth', 1.5);

% 选取 N=64, 80, 96, 112, 128 绘制轨迹（若 nN<25 则自适应选取）
idx_plot = round(linspace(1, nN, min(5, nN)));
colors = lines(numel(idx_plot));
leg_entries = {'BS', 'Users', 'User Region', '初始轨迹'};
for j = 1:numel(idx_plot)
    ii = idx_plot(j);
    plot(x_opt_cell{ii}, y_opt_cell{ii}, '-', 'LineWidth', 1.5, 'Color', colors(j,:));
    leg_entries{end+1} = sprintf('N=%d', N_list(ii));
end
legend(leg_entries, 'Location', 'best');
title('有源 RIS：不同 N 下的优化轨迹');
xlabel('x (m)'); ylabel('y (m)');

% 保存结果
save('study_N_results.mat', 'N_list', 'Rsum_list', 'wsr_mean_list', 'wsr_max_list', 'x_opt_cell', 'y_opt_cell', 'rate_hist_cell', 'user_pos', 'bs_pos', 'x_traj', 'y_traj');
fprintf('\n========== 研究完成，结果已保存至 study_N_results.mat ==========\n');
