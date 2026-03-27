clear; clc; close all;

% =========================
% 实验开关
% =========================
E_test_mode = true;
RIS_ablation_mode = true;

% A组：不改直射链路强度，只画出逐slot的 Active - NoRIS 增益
A_gain_curve_mode = true;

% B组：引入“每个slot一个小尺度衰落 realization”
% 同一slot内部 inner-iteration 固定，不同slot不同
B_slot_fading_mode = true;

% =========================
% Scenario setup
% =========================
x_start = 500;
y_start = -300;
y_min = -300;
y_max = 100;

v_max = 10;
T_total = 480;
Nslot = 100;
dt = T_total / Nslot;

if E_test_mode
    Nslot = 160;   % shorten runtime for testing
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
P_total_mW = 10^(P_total_dBm/10);
RIS_power_ratio = 0.30;
Pr_max = P_total_mW * RIS_power_ratio;
Ps_max = P_total_mW * (1 - RIS_power_ratio);

sigma2 = 1e-13;
sigmar2 = 1e-13;
eta_k = ones(K,1);
f_c = 5;
lambda = 3e8/(f_c*1e9);
K_rician = 10;
tau_power = 4;
kappa_power = 4;
gamma_reflect = 0;
scale_G = 1e3;
scale_f = 1e6;

x_traj = x_start * ones(Nslot, 1);
y_traj = y_start * ones(Nslot, 1);

% =========================
% Geometry
% =========================
bs_pos = [0, 0, ht_BS];

rng(42);
theta = 2*pi*rand(K,1);
rad = user_R*rand(K,1);
user_pos = zeros(K,3);
for k = 1:K
    user_pos(k,:) = [ ...
        user_X_center + rad(k)*cos(theta(k)), ...
        user_Y_center + rad(k)*sin(theta(k)), ...
        hr_User];
end

% =========================
% Initial communication variables
% =========================
rng(123);
h_k0 = zeros(M,K);
f_k0 = zeros(N,K);
G0 = zeros(N,M);
Theta0 = diag(exp(1j*2*pi*rand(N,1)));
W0 = exp(1j*2*pi*rand(K*M,1))*sqrt(Ps_max/K/M);

% =========================
% Trajectory config
% =========================
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

% ===== B组新增：每个slot一个固定seed =====
traj_cfg.use_slot_seed = B_slot_fading_mode;
traj_cfg.slot_seed_base = 1000;

fprintf('========== Slot-by-slot USV optimization ==========\n');
fprintf('A_gain_curve_mode    = %d\n', A_gain_curve_mode);
fprintf('B_slot_fading_mode   = %d\n', B_slot_fading_mode);

% =========================
% Run slot-by-slot optimization
% =========================
if E_test_mode
    cfg_B = traj_cfg;
    cfg_B.e2_ratio_cap = 5.0;
    cfg_B.e_clip = 10.0;
    [W_opt, Theta_opt, Rsum, ~, ~, x_opt, y_opt, rate_hist, ~] = active_RIS_precoding( ...
        M, K, N, Ps_max, Pr_max, sigma2, sigmar2, eta_k, ...
        Theta0, W0, h_k0, f_k0, G0, cfg_B);
else
    [W_opt, Theta_opt, Rsum, ~, ~, x_opt, y_opt, rate_hist, ~] = active_RIS_precoding( ...
        M, K, N, Ps_max, Pr_max, sigma2, sigmar2, eta_k, ...
        Theta0, W0, h_k0, f_k0, G0, traj_cfg);
end

wsr_active = rate_hist(:);
[wsr_max, best_slot] = max(wsr_active);

% =========================
% A组：逐slot NoRIS基线 + 增益曲线
% =========================
wsr_noRIS = nan(Nslot, 1);
gain_abs  = nan(Nslot, 1);
gain_pct  = nan(Nslot, 1);

sigma2_eff  = sigma2 * (scale_G^2 * scale_f^2);
sigmar2_eff = sigmar2 * (scale_G^2);
Pr_max_eff  = Pr_max * (scale_G^2);

if A_gain_curve_mode
    fprintf('\n========== A组：逐slot Active vs NoRIS 增益评估 ==========\n');
    for n = 1:Nslot
        if B_slot_fading_mode
            slot_seed = traj_cfg.slot_seed_base + n;
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

        h_n = 1 * ch_n(1).h_k;  % 与主流程保持一致
        f_n = ch_n(1).f_k;
        G_n = ch_n(1).G;

        [~, wsr_noRIS(n)] = NoRIS_precoding( ...
            M, K, N, Ps_max, sigma2_eff, eta_k, W0, h_n, f_n, G_n);

        gain_abs(n) = wsr_active(n) - wsr_noRIS(n);
        gain_pct(n) = 100 * gain_abs(n) / max(wsr_noRIS(n), 1e-9);

        fprintf('slot=%3d | Active=%.4f | NoRIS=%.4f | Gain=%.4f (%.2f%%)\n', ...
            n, wsr_active(n), wsr_noRIS(n), gain_abs(n), gain_pct(n));
    end
end

% =========================
% 保留 best-slot 消融输出
% =========================
if RIS_ablation_mode
    if A_gain_curve_mode
        fprintf('\nRIS ablation @best slot %d: ActiveRIS=%.4f, NoRIS=%.4f, Gain=%.4f (%.2f%%)\n', ...
            best_slot, wsr_active(best_slot), wsr_noRIS(best_slot), ...
            gain_abs(best_slot), gain_pct(best_slot));
    else
        if B_slot_fading_mode
            slot_seed = traj_cfg.slot_seed_base + best_slot;
            [ch_best, ~, ~] = build_slot_channels( ...
                K, N, M, bs_pos, user_pos, x_opt(best_slot), y_opt(best_slot), ...
                h_RIS, lambda, gamma_reflect, K_rician, false, ...
                scale_G, scale_f, slot_seed);
        else
            [ch_best, ~, ~] = build_slot_channels( ...
                K, N, M, bs_pos, user_pos, x_opt(best_slot), y_opt(best_slot), ...
                h_RIS, lambda, gamma_reflect, K_rician, false, ...
                scale_G, scale_f);
        end

        h_best = 1 * ch_best(1).h_k;
        f_best = ch_best(1).f_k;
        G_best = ch_best(1).G;

        [~, ~, wsr_active_best, ~, ~] = active_RIS_precoding( ...
            M, K, N, Ps_max, Pr_max_eff, sigma2_eff, sigmar2_eff, ...
            eta_k, Theta0, W0, h_best, f_best, G_best);

        [~, wsr_noRIS_best] = NoRIS_precoding( ...
            M, K, N, Ps_max, sigma2_eff, eta_k, W0, h_best, f_best, G_best);

        gain_abs_best = wsr_active_best - wsr_noRIS_best;
        gain_pct_best = 100 * gain_abs_best / max(wsr_noRIS_best, 1e-9);

        fprintf('\nRIS ablation @best slot %d: ActiveRIS=%.4f, NoRIS=%.4f, Gain=%.4f (%.2f%%)\n', ...
            best_slot, wsr_active_best, wsr_noRIS_best, gain_abs_best, gain_pct_best);
    end
end

% =========================
% Plot
% =========================
figure('Position', [80, 80, 1350, 760]);

% ---- subplot 1: trajectory ----
subplot(2, 2, 1); hold on; grid on; box on;
plot(bs_pos(1), bs_pos(2), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
scatter(user_pos(:,1), user_pos(:,2), 50, 'b', 'filled');

th = linspace(0, 2*pi, 100);
plot(user_X_center + user_R*cos(th), user_Y_center + user_R*sin(th), '--', 'Color', [0.6 0.6 0.6]);
plot(x_opt, y_opt, '-r', 'LineWidth', 2, 'Marker', '.');

step_mark = 5;
for n = 1:step_mark:Nslot
    y_offset = 20;
    text(x_opt(n), y_opt(n) + y_offset, sprintf('t=%d', n), ...
        'FontSize', 9, 'Color', [0.5 0 0.5]);
    plot(x_opt(n), y_opt(n), 'o', 'MarkerSize', 5, ...
        'MarkerEdgeColor', [0.5 0 0.5]);
end

legend({'BS', 'Users', 'User Region', 'Active RIS Trajectory'}, 'Location', 'best');
title('USV Trajectory (Slot-by-slot)');
xlabel('x (m)');
ylabel('y (m)');

% ---- subplot 2: active vs noRIS ----
subplot(2, 2, 2); hold on; grid on; box on;
plot(1:Nslot, wsr_active, '-r', 'LineWidth', 2.0, 'Marker', 'o', ...
    'MarkerSize', 3, 'MarkerFaceColor', 'r');

if A_gain_curve_mode
    plot(1:Nslot, wsr_noRIS, '-k', 'LineWidth', 1.8, 'Marker', 's', ...
        'MarkerSize', 3);
    legend({'Active RIS', 'No RIS'}, 'Location', 'best');
else
    legend({'Active RIS'}, 'Location', 'best');
end

plot(best_slot, wsr_max, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
text(best_slot, wsr_max + 0.5, sprintf('Active Max: %.2f', wsr_max), ...
    'HorizontalAlignment', 'center');

xlabel('Slot Index');
ylabel('WSR (bps/Hz)');
title(sprintf('WSR Curve (Active Max: %.4f)', wsr_max));

% ---- subplot 3: absolute gain ----
subplot(2, 2, 3); hold on; grid on; box on;
if A_gain_curve_mode
    [gain_max_abs, gain_best_slot] = max(gain_abs);
    plot(1:Nslot, gain_abs, '-b', 'LineWidth', 2.0, 'Marker', 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', 'b');
    plot(gain_best_slot, gain_max_abs, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
    text(gain_best_slot, gain_max_abs + 0.2, sprintf('Max Gain: %.2f', gain_max_abs), ...
        'HorizontalAlignment', 'center');
end
xlabel('Slot Index');
ylabel('\DeltaWSR = Active - NoRIS');
title('A组：Absolute Gain of Active RIS');

% ---- subplot 4: percentage gain ----
subplot(2, 2, 4); hold on; grid on; box on;
if A_gain_curve_mode
    [gain_max_pct, gain_best_slot_pct] = max(gain_pct);
    plot(1:Nslot, gain_pct, '-m', 'LineWidth', 2.0, 'Marker', 'o', ...
        'MarkerSize', 3, 'MarkerFaceColor', 'm');
    plot(gain_best_slot_pct, gain_max_pct, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
    text(gain_best_slot_pct, gain_max_pct + 0.5, sprintf('Max Gain: %.2f%%', gain_max_pct), ...
        'HorizontalAlignment', 'center');
end
xlabel('Slot Index');
ylabel('Gain (%)');
title('A组：Relative Gain of Active RIS');

sgtitle('Slot-by-slot USV Optimization + A组增益分析 + B组逐slot小尺度起伏');