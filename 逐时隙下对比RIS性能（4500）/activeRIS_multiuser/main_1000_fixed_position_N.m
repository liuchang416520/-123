clear; clc; close all;
% ==========================================================
% 有源 RIS 固定位置 N 敏感性研究
% 1) 先做轨迹优化，找到 WSR 最高的 RIS 位置
% 2) 固定该位置，研究反射元件数量 N 对 WSR 的影响
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

% 5. 高度与通信设置
ht_BS = 5; h_RIS = 2; hr_User = 1.5;
K = 4; M = 10; N_base = 64;
P_total_dBm = 15;
P_total_mW  = 10^(P_total_dBm/10);
RIS_power_ratio = 0.1;
Pr_max = P_total_mW * RIS_power_ratio;
Ps_max = P_total_mW * (1 - RIS_power_ratio);

sigma2 = 1e-13; sigmar2 = 1e-13; eta_k = ones(K,1);
f_c = 5; lambda = 3e8/(f_c*1e9); K_rician = 10;
tau_power = 4; kappa_power = 4; gamma_reflect = 0;
scale_G = 1e3; scale_f = 1e6;

% 有效功率/噪声（与轨迹优化一致）
sigma2_eff = sigma2 * (scale_G^2 * scale_f^2);
sigmar2_eff = sigmar2 * (scale_G^2);
Pr_max_eff = Pr_max * (scale_G^2);

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

% ==========================================================
% 阶段一：有源 RIS 轨迹优化，找到 WSR 最高的位置
% ==========================================================
fprintf('========== 阶段一：有源 RIS 轨迹优化 (N=%d) ==========\n', N_base);

rng(123);
h_k0 = zeros(M, K); f_k0 = zeros(N_base, K); G0 = zeros(N_base, M);
Theta0 = diag(exp(1j*2*pi*rand(N_base,1)));
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
traj_cfg.no_plot = true;  % 由本脚本统一绘图

[~, ~, ~, ~, ~, x_opt, y_opt, rate_hist, ~] = active_RIS_precoding(M, K, N_base, Ps_max, Pr_max, sigma2, sigmar2, eta_k, Theta0, W0, h_k0, f_k0, G0, traj_cfg);

wsr_final = rate_hist(:, end);
[wsr_max_val, best_slot] = max(wsr_final);
x_best = x_opt(best_slot);
y_best = y_opt(best_slot);

fprintf('  WSR 最高时隙: %d, 位置 (%.1f, %.1f) m, WSR=%.4f bps/Hz\n', best_slot, x_best, y_best, wsr_max_val);

% ==========================================================
% 阶段二：固定该位置，扫描 N=64:2:128，研究 WSR 随 N 变化
% ==========================================================
fprintf('\n========== 阶段二：固定位置 (%.1f, %.1f)，N 敏感性研究 ==========\n', x_best, y_best);

N_list = 64:100:864;
nN = numel(N_list);
Rsum_fixed = zeros(nN, 1);

for ii = 1:nN
    N = N_list(ii);
    fprintf('  N = %d (%d/%d) ... ', N, ii, nN);

    % 在固定位置构建信道（维度随 N 变化）
    [channels, ~, ~] = build_slot_channels(K, N, M, bs_pos, user_pos, [x_best], [y_best], h_RIS, lambda, gamma_reflect, K_rician, false, scale_G, scale_f);
    ch = channels(1);
    ch.h_k = 0.3 * ch.h_k;  % 与轨迹优化一致的直射衰减

    % 初值
    rng(123 + ii);  % 不同 N 用不同初值，避免相关性
    Theta0_n = diag(exp(1j*2*pi*rand(N,1)));
    W0_n = exp(1j*2*pi*rand(K*M,1)) * sqrt(Ps_max/K/M);

    % 单点通信优化（无轨迹）
    [~, ~, Rsum_fixed(ii)] = active_RIS_precoding(M, K, N, Ps_max, Pr_max_eff, sigma2_eff, sigmar2_eff, eta_k, Theta0_n, W0_n, ch.h_k, ch.f_k, ch.G);

    fprintf('WSR=%.4f\n', Rsum_fixed(ii));
end

% ==========================================================
% 结果绘图
% ==========================================================
figure('Position', [50, 50, 1100, 450]);

% 子图一：轨迹图，标出 WSR 最高位置
subplot(1, 2, 1); hold on; grid on; box on;
plot(bs_pos(1), bs_pos(2), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
scatter(user_pos(:,1), user_pos(:,2), 50, 'b', 'filled');
th = linspace(0, 2*pi, 100);
plot(user_X_center + user_R*cos(th), user_Y_center + user_R*sin(th), '--', 'Color', [0.6 0.6 0.6]);
plot(x_traj, y_traj, '--k', 'LineWidth', 1.5);
plot(x_opt, y_opt, '-b', 'LineWidth', 1.5);
plot(x_best, y_best, 'rp', 'MarkerSize', 16, 'MarkerFaceColor', 'r', 'LineWidth', 2);
text(x_best, y_best + 25, sprintf('WSR最高\n(%.0f,%.0f)', x_best, y_best), 'HorizontalAlignment', 'center', 'FontSize', 9);
legend({'BS', 'Users', 'User Region', '初始轨迹', '优化轨迹', '固定研究位置'}, 'Location', 'best');
title(sprintf('有源 RIS 轨迹 (WSR最高: 时隙%d, (%.0f,%.0f)m)', best_slot, x_best, y_best));
xlabel('x (m)'); ylabel('y (m)');

% 子图二：固定位置下 WSR vs N
subplot(1, 2, 2); hold on; grid on; box on;
plot(N_list, Rsum_fixed, '-b', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
xlabel('反射元件数量 N');
ylabel('WSR (bps/Hz)');
title(sprintf('固定位置 (%.0f, %.0f) m：WSR 随 N 变化', x_best, y_best));

% 保存结果
save('fixed_position_N_results.mat', 'x_best', 'y_best', 'best_slot', 'N_list', 'Rsum_fixed', 'x_opt', 'y_opt', 'wsr_final');
fprintf('\n========== 研究完成，结果已保存至 fixed_position_N_results.mat ==========\n');
