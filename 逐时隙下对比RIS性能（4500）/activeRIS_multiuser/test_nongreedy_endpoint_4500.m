clear; clc; close all;
% ==========================================================
% Test script: non-greedy slot-by-slot AO/SCA with hard terminal point
% Final point is forced to (4500, -300).
% ==========================================================

% 1) start / terminal points
x_start = 500;     y_start = -300;
x_end   = 4500;   y_end   = -300;

% 2) boundaries and motion
v_max = 10;
T_total = 800;
Nslot = 100;
dt = T_total / Nslot;
y_min = -300;
y_max = 0;

% 3) users (static)
user_X_center = 4500;
user_Y_center = 0;
user_R = 100;

% 4) heights and comm params
ht_BS = 10; h_RIS = 1.5; hr_User = 3;
K = 4; M = 10; N = 64;
P_total_dBm = 15;
P_total_mW  = 10^(P_total_dBm/10);
RIS_power_ratio = 0.1;
Pr_max = P_total_mW * RIS_power_ratio;
Ps_max = P_total_mW * (1 - RIS_power_ratio);

sigma2 = 1e-13; sigmar2 = 1e-13; eta_k = ones(K,1);
f_c = 5; lambda = 3e8/(f_c*1e9); K_rician = 10;
tau_power = 4; kappa_power = 4; gamma_reflect = 0;
scale_G = 1e3; scale_f = 1e6;

% 5) initial trajectory: straight line start -> end
x_traj = linspace(x_start, x_end, Nslot).';
y_traj = linspace(y_start, y_end, Nslot).';

bs_pos = [0, 0, ht_BS];
rng(42);
theta_u = 2*pi*rand(K,1); rad = user_R*rand(K,1);
user_pos = zeros(K,3);
for k = 1:K
    user_pos(k,:) = [user_X_center + rad(k)*cos(theta_u(k)), user_Y_center + rad(k)*sin(theta_u(k)), hr_User];
end

% 6) init channels / vars
rng(123);
h_k0 = zeros(M,K); f_k0 = zeros(N,K); G0 = zeros(N,M);
Theta0 = diag(exp(1j*2*pi*rand(N,1)));
W0 = exp(1j*2*pi*rand(K*M,1))*sqrt(Ps_max/K/M);

% 7) trajectory config
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
traj_cfg.lambda = lambda;
traj_cfg.gamma_reflect = gamma_reflect;
traj_cfg.K_rician = K_rician;
traj_cfg.tau_power = tau_power;
traj_cfg.kappa_power = kappa_power;
traj_cfg.no_direct_link = false;
traj_cfg.MaxIter = 12;
traj_cfg.MinIter = 3;
traj_cfg.eps_slot = 1e-3;
traj_cfg.scale_G = scale_G;
traj_cfg.scale_f = scale_f;
traj_cfg.use_smart_init = false;
traj_cfg.trust_step = 0.8 * v_max * dt;
traj_cfg.log_level = 1;

% strict terminal design
traj_cfg.q_final = [x_end; y_end];
traj_cfg.fix_endpoint_last_slot = true;
traj_cfg.endpoint_reachability = true;
traj_cfg.lambda_end = 5e-5;   % soft pull before last slot
traj_cfg.L_sca = 1.0;
traj_cfg.x_coast_min = 0;
traj_cfg.safe_user_r = 20 * ones(K,1);   % optional

[W_opt, Theta_opt, Rsum, ~, ~, x_opt, y_opt, rate_hist, slot_info] = ...
    active_RIS_precoding_nongreedy(M, K, N, Ps_max, Pr_max, sigma2, sigmar2, eta_k, ...
    Theta0, W0, h_k0, f_k0, G0, traj_cfg);

fprintf('\nFinal position = (%.2f, %.2f)\n', x_opt(end), y_opt(end));
fprintf('Required final  = (%.2f, %.2f)\n', x_end, y_end);
fprintf('Final slot WSR  = %.4f\n', rate_hist(end));

figure('Position',[100,100,1100,420]);
subplot(1,2,1); hold on; grid on; box on;
plot(bs_pos(1), bs_pos(2), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
scatter(user_pos(:,1), user_pos(:,2), 50, 'b', 'filled');
plot(x_traj, y_traj, '--k', 'LineWidth', 1.2);
plot(x_opt, y_opt, '-r', 'LineWidth', 2);
plot(x_start, y_start, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot(x_end, y_end, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
legend({'BS','Users','Init','Optimized','Start','End'}, 'Location', 'best');
xlabel('x (m)'); ylabel('y (m)'); title('Non-greedy slot-by-slot trajectory');

subplot(1,2,2); hold on; grid on; box on;
plot(1:Nslot, rate_hist, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('slot'); ylabel('WSR (bps/Hz)'); title('Slot WSR');
