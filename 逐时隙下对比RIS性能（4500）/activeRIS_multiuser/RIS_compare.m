clear; clc; close all;

% ==========================================================
% main_1000.m
% 目标：
% 1) 先用 Active RIS 的逐时隙轨迹优化得到公共轨迹 x_opt / y_opt
% 2) 在这条公共轨迹上，对每个 slot 的同一组信道分别计算：
%       - Active RIS
%       - Passive RIS
%       - No RIS
% 3) 输出三条 WSR 曲线，以及：
%       - Active - NoRIS
%       - Active - Passive
%
% 重要前提（请先确认）：
% A. build_slot_channels.m 里已经统一设置了直射链路缩放
% B. active_RIS_precoding.m 里额外那两处 h_k 手动缩放已经改成 1 或删除
% 否则 Active / Passive / NoRIS 看到的直射链路不一致，比较不公平
% ==========================================================

%% =========================
% 开关
% ==========================
E_test_mode = true;          % true: Nslot=160 便于快速测试
B_slot_fading_mode = true;   % 若你已按前面建议改过 build_slot_channels，可开启逐slot realization
slot_seed_base = 1000;

% 这里假设直射链路已经在 build_slot_channels.m 统一缩放过
% 因此这里不再做第二次缩放，统一乘 1
DIRECT_LINK_POST_SCALE = 1.0;

%% =========================
% 场景参数
% ==========================
x_start = 500;
y_start = -300;
y_min = -300;
y_max = 100;

v_max = 10;
T_total = 480;
Nslot = 100;
dt = T_total / Nslot;

if E_test_mode
    Nslot = 160;
    dt = T_total / Nslot;
end

user_X_center = 4500;
user_Y_center = 0;
user_R = 100;

ht_BS = 10;
h_RIS = 3;
hr_User = 3;

K = 4;
M = 10;
N = 128;

P_total_dBm = 20;
P_total_mW = 10^(P_total_dBm / 10);
RIS_power_ratio = 0.10;

% --------------------------
% 功率公平规则（沿用 main_fixed_RIS_812_N_study.m）
% Active：BS 功率减少，RIS 单独拿功率
% Passive / NoRIS：BS 独占总功率
% --------------------------
Pr_max = P_total_mW * RIS_power_ratio;
Ps_max_active = P_total_mW * (1 - RIS_power_ratio);
Ps_max_passive = P_total_mW;
Ps_max_noRIS = P_total_mW;

sigma2 = 1e-13;
sigmar2 = 1e-13;
eta_k = ones(K,1);

f_c = 5;
lambda = 3e8 / (f_c * 1e9);
K_rician = 10;

tau_power = 4;
kappa_power = 4;
gamma_reflect = 0;

scale_G = 1e3;
scale_f = 1e6;

sigma2_eff  = sigma2  * (scale_G^2 * scale_f^2);
sigmar2_eff = sigmar2 * (scale_G^2);
Pr_max_eff  = Pr_max  * (scale_G^2);

%% =========================
% 初始轨迹
% ==========================
x_traj = x_start * ones(Nslot, 1);
y_traj = y_start * ones(Nslot, 1);

%% =========================
% 几何场景
% ==========================
bs_pos = [0, 0, ht_BS];

rng(42);
theta_user = 2*pi*rand(K,1);
rad_user = user_R * rand(K,1);
user_pos = zeros(K,3);
for k = 1:K
    user_pos(k,:) = [ ...
        user_X_center + rad_user(k)*cos(theta_user(k)), ...
        user_Y_center + rad_user(k)*sin(theta_user(k)), ...
        hr_User];
end

%% =========================
% Active 轨迹优化的初始化
% ==========================
rng(123);
h_k0 = zeros(M, K);
f_k0 = zeros(N, K);
G0   = zeros(N, M);

Theta0_active_traj = diag(exp(1j * 2*pi * rand(N,1)));
W0_active_traj = exp(1j * 2*pi * rand(K*M,1)) * sqrt(Ps_max_active / K / M);

%% =========================
% Active 轨迹优化配置
% ==========================
traj_cfg = struct();
traj_cfg.enable = true;
traj_cfg.bs_pos = bs_pos;
traj_cfg.user_pos = user_pos;
traj_cfg.x_traj = x_traj;
traj_cfg.y_traj = y_traj;
traj_cfg.h_const = h_RIS;
traj_cfg.Nslot = Nslot;
traj_cfg.dt = dt;
traj_cfg.v_max = v_max;
traj_cfg.y_min = y_min;
traj_cfg.y_max = y_max;
traj_cfg.x_coast_min = 0;
traj_cfg.lambda = lambda;
traj_cfg.gamma_reflect = gamma_reflect;
traj_cfg.K_rician = K_rician;
traj_cfg.tau_power = tau_power;
traj_cfg.kappa_power = kappa_power;
traj_cfg.no_direct_link = false;
traj_cfg.MaxIter = 16;
traj_cfg.MinIter = 5;
traj_cfg.disp_thresh = 15;
traj_cfg.return_to_endpoint = false;
traj_cfg.enable_early_hover_stop = false;
traj_cfg.scale_G = scale_G;
traj_cfg.scale_f = scale_f;
traj_cfg.use_smart_init = true;
traj_cfg.disable_cross_slot_warm_start_when_no_direct = true;
traj_cfg.zero_wsr_thresh = 1e-6;
traj_cfg.verbose_comm = true;
traj_cfg.log_level = 1;

% 若你已经按前面建议改过 build_slot_channels / active_RIS_precoding，可开启
traj_cfg.use_slot_seed = B_slot_fading_mode;
traj_cfg.slot_seed_base = slot_seed_base;

% E 平衡参数
traj_cfg.e2_ratio_cap = 5.0;
traj_cfg.e_clip = 10.0;

fprintf('========== Step 1: Active RIS slot-by-slot trajectory optimization ==========\n');

[~, ~, ~, ~, ~, x_opt, y_opt, rate_hist_active_traj, ~] = active_RIS_precoding( ...
    M, K, N, Ps_max_active, Pr_max, sigma2, sigmar2, eta_k, ...
    Theta0_active_traj, W0_active_traj, h_k0, f_k0, G0, traj_cfg);

%% =========================
% Step 2: 在同一条公共轨迹上做三路通信对比
% ==========================
fprintf('\n========== Step 2: Same-trajectory communication comparison ==========\n');

% 三条曲线
wsr_active  = zeros(Nslot, 1);
wsr_passive = zeros(Nslot, 1);
wsr_noRIS   = zeros(Nslot, 1);

% 增益
gain_active_over_noRIS   = zeros(Nslot, 1);
gain_active_over_passive = zeros(Nslot, 1);
gain_pct_over_noRIS      = zeros(Nslot, 1);
gain_pct_over_passive    = zeros(Nslot, 1);

% 为三条路分别维护独立的热启动状态，避免互相污染
rng(2026);
Theta_active_eval = diag(exp(1j * 2*pi * rand(N,1)));
W_active_eval = exp(1j * 2*pi * rand(K*M,1)) * sqrt(Ps_max_active / K / M);

rng(2027);
Theta_passive_eval = diag(exp(1j * 2*pi * rand(N,1)));
W_passive_eval = exp(1j * 2*pi * rand(K*M,1)) * sqrt(Ps_max_passive / K / M);

rng(2028);
W_noRIS_eval = exp(1j * 2*pi * rand(K*M,1)) * sqrt(Ps_max_noRIS / K / M);

for n = 1:Nslot
    % --------------------------
    % 同一slot只生成一次公共信道
    % --------------------------
    if B_slot_fading_mode
        slot_seed = slot_seed_base + n;
        [ch_n, ~, ~] = build_slot_channels( ...
            K, N, M, bs_pos, user_pos, x_opt(n), y_opt(n), ...
            h_RIS, lambda, gamma_reflect, K_rician, false, ...
            scale_G, scale_f, slot_seed);
    else
        [ch_n, ~, ~] = build_slot_channels( ...
            K, N, M, bs_pos, user_pos, x_opt(n), y_opt(n), ...
            h_RIS, lambda, gamma_reflect, K_rician, false, ...
            scale_G, scale_f);
    end

    % 这里不再做额外不对称缩放
    h_n = DIRECT_LINK_POST_SCALE * ch_n(1).h_k;
    f_n = ch_n(1).f_k;
    G_n = ch_n(1).G;

    % --------------------------
    % No RIS
    % --------------------------
    [W_noRIS_eval, wsr_noRIS(n)] = NoRIS_precoding( ...
        M, K, N, Ps_max_noRIS, sigma2_eff, eta_k, ...
        W_noRIS_eval, h_n, f_n, G_n);

    % --------------------------
    % Passive RIS
    % --------------------------
    [W_passive_eval, Theta_passive_eval, wsr_passive(n)] = passive_RIS_precoding( ...
        M, K, N, Ps_max_passive, sigma2_eff, eta_k, ...
        Theta_passive_eval, W_passive_eval, h_n, f_n, G_n);

    % --------------------------
    % Active RIS
    % 这里用固定位置通信优化模式重新评估，
    % 以便和 Passive / NoRIS 一样都在同一公共信道上比较
    % --------------------------
    [W_active_eval, Theta_active_eval, wsr_active(n), ~, ~] = active_RIS_precoding( ...
        M, K, N, Ps_max_active, Pr_max_eff, sigma2_eff, sigmar2_eff, eta_k, ...
        Theta_active_eval, W_active_eval, h_n, f_n, G_n);

    gain_active_over_noRIS(n)   = wsr_active(n) - wsr_noRIS(n);
    gain_active_over_passive(n) = wsr_active(n) - wsr_passive(n);

    gain_pct_over_noRIS(n) = 100 * gain_active_over_noRIS(n) / max(wsr_noRIS(n), 1e-9);
    gain_pct_over_passive(n) = 100 * gain_active_over_passive(n) / max(wsr_passive(n), 1e-9);

    fprintf(['slot=%3d | Active=%.4f | Passive=%.4f | NoRIS=%.4f | ' ...
             'A-P=%.4f | A-N=%.4f\n'], ...
        n, wsr_active(n), wsr_passive(n), wsr_noRIS(n), ...
        gain_active_over_passive(n), gain_active_over_noRIS(n));
end

[wsr_active_max, best_slot_active] = max(wsr_active);
[gap_noRIS_max, best_slot_gap_noRIS] = max(gain_active_over_noRIS);
[gap_passive_max, best_slot_gap_passive] = max(gain_active_over_passive);

%% =========================
% 绘图
% ==========================
figure('Position', [80, 80, 1400, 760]);

% ---- 1. 轨迹图 ----
subplot(2, 2, 1); hold on; grid on; box on;
plot(bs_pos(1), bs_pos(2), 'ks', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'g', 'LineWidth', 1.5);
scatter(user_pos(:,1), user_pos(:,2), 50, 'b', 'filled');

th = linspace(0, 2*pi, 100);
plot(user_X_center + user_R*cos(th), ...
     user_Y_center + user_R*sin(th), '--', 'Color', [0.6 0.6 0.6]);

plot(x_opt, y_opt, '-r', 'LineWidth', 2, 'Marker', '.');

step_mark = 5;
for n = 1:step_mark:Nslot
    text(x_opt(n), y_opt(n) + 20, sprintf('t=%d', n), ...
        'FontSize', 8, 'Color', [0.5 0 0.5]);
    plot(x_opt(n), y_opt(n), 'o', 'MarkerSize', 4, ...
        'MarkerEdgeColor', [0.5 0 0.5]);
end

legend({'BS', 'Users', 'User Region', 'Public Trajectory'}, 'Location', 'best');
title('Common Trajectory (optimized by Active RIS)');
xlabel('x (m)');
ylabel('y (m)');

% ---- 2. 三路 WSR 对比 ----
subplot(2, 2, 2); hold on; grid on; box on;
plot(1:Nslot, wsr_active,  '-r', 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 3);
plot(1:Nslot, wsr_passive, '-b', 'LineWidth', 2.0, 'Marker', 's', 'MarkerSize', 3);
plot(1:Nslot, wsr_noRIS,   '-k', 'LineWidth', 2.0, 'Marker', '^', 'MarkerSize', 3);

plot(best_slot_active, wsr_active_max, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
text(best_slot_active, wsr_active_max + 0.5, sprintf('Active Max: %.2f', wsr_active_max), ...
    'HorizontalAlignment', 'center');

xlabel('Slot Index');
ylabel('WSR (bps/Hz)');
title('Active vs Passive vs NoRIS on the Same Trajectory');
legend({'Active RIS', 'Passive RIS', 'No RIS', 'Active Max'}, 'Location', 'best');

% ---- 3. Active 相对 NoRIS 的增益 ----
subplot(2, 2, 3); hold on; grid on; box on;
yyaxis left;
plot(1:Nslot, gain_active_over_noRIS, '-m', 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 3);
ylabel('\DeltaWSR = Active - NoRIS');

yyaxis right;
plot(1:Nslot, gain_pct_over_noRIS, '--c', 'LineWidth', 1.8);
ylabel('Gain over NoRIS (%)');

plot(best_slot_gap_noRIS, gain_active_over_noRIS(best_slot_gap_noRIS), ...
    'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');

xlabel('Slot Index');
title('Active over NoRIS');

% ---- 4. Active 相对 Passive 的增益 ----
subplot(2, 2, 4); hold on; grid on; box on;
yyaxis left;
plot(1:Nslot, gain_active_over_passive, '-g', 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 3);
ylabel('\DeltaWSR = Active - Passive');

yyaxis right;
plot(1:Nslot, gain_pct_over_passive, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8);
ylabel('Gain over Passive (%)');

plot(best_slot_gap_passive, gain_active_over_passive(best_slot_gap_passive), ...
    'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');

xlabel('Slot Index');
title('Active over Passive');

sgtitle('Same-condition Comparison: Active / Passive / NoRIS');

%% =========================
% 控制台摘要
% ==========================
fprintf('\n========== Summary ==========\n');
fprintf('Active   : mean=%.4f, max=%.4f @ slot %d\n', mean(wsr_active), max(wsr_active), best_slot_active);
fprintf('Passive  : mean=%.4f, max=%.4f\n', mean(wsr_passive), max(wsr_passive));
fprintf('NoRIS    : mean=%.4f, max=%.4f\n', mean(wsr_noRIS), max(wsr_noRIS));

fprintf('Active - NoRIS   : mean=%.4f, max=%.4f @ slot %d, mean%%=%.2f%%\n', ...
    mean(gain_active_over_noRIS), max(gain_active_over_noRIS), ...
    best_slot_gap_noRIS, mean(gain_pct_over_noRIS));

fprintf('Active - Passive : mean=%.4f, max=%.4f @ slot %d, mean%%=%.2f%%\n', ...
    mean(gain_active_over_passive), max(gain_active_over_passive), ...
    best_slot_gap_passive, mean(gain_pct_over_passive));

save('main_1000_three_way_compare.mat', ...
    'x_opt', 'y_opt', ...
    'wsr_active', 'wsr_passive', 'wsr_noRIS', ...
    'gain_active_over_noRIS', 'gain_active_over_passive', ...
    'gain_pct_over_noRIS', 'gain_pct_over_passive', ...
    'rate_hist_active_traj');

fprintf('Results saved to main_1000_three_way_compare.mat\n');