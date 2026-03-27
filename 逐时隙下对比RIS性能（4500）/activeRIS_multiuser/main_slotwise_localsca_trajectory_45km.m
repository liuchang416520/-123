clear; clc; close all;

%% =========================================
% 1000m 场景：逐时隙 AO + 局部 SCA 主脚本
% 目标：在每个时隙内交替优化通信变量与位置
% 并通过回溯线搜索保证“WSR 不下降”
% ==========================================

%% ========== 场景参数 ==========
x_start = 3000;
y_start = -300;
y_min = -300;
y_max = 100;

% x 方向活动范围（避免 USV 无约束漂移）
x_min = 0;
x_max = 6000;

v_max = 10;       % m/s
T_total = 250;    % s
Nslot = 100;
dt = T_total / Nslot;

user_X_center = 4500;
user_Y_center = 0;
user_R = 100;

ht_BS = 10;
h_RIS = 3;
hr_User = 3;

K = 4;      % 用户数
M = 10;     % BS 天线数
N = 128;    % RIS 单元数

P_total_dBm = 25;
P_total_mW = 10^(P_total_dBm/10);
RIS_power_ratio = 0.1;
Pr_max = P_total_mW * RIS_power_ratio;
Ps_max = P_total_mW * (1 - RIS_power_ratio);

sigma2  = 1e-13;
sigmar2 = 1e-13;
eta_k   = ones(K,1);

f_c = 5;  % GHz
lambda = 3e8/(f_c*1e9);
K_rician = 10;
gamma_reflect = 0;

% 与主流程保持一致的缩放
scale_G = 1e3;
scale_f = 1e6;

%% ========== 几何初始化 ==========
bs_pos = [0, 0, ht_BS];

rng(42);
theta = 2*pi*rand(K,1);
rad   = user_R*rand(K,1);
user_pos = zeros(K,3);
for k = 1:K
    user_pos(k,:) = [user_X_center + rad(k)*cos(theta(k)), ...
                     user_Y_center + rad(k)*sin(theta(k)), ...
                     hr_User];
end

fprintf('========== 1000m 场景：逐时隙 AO + 局部 SCA ==========\n');
fprintf('BS 位置 = [%.2f, %.2f, %.2f]\n', bs_pos(1), bs_pos(2), bs_pos(3));
for k = 1:K
    fprintf('User %d 位置 = [%.2f, %.2f, %.2f]\n', ...
        k, user_pos(k,1), user_pos(k,2), user_pos(k,3));
end

%% ========== 通信变量初始化 ==========
rng(123);
Theta_cur = diag(exp(1j*2*pi*rand(N,1)));
W_cur = exp(1j*2*pi*rand(K*M,1))*sqrt(Ps_max/K/M);

% 占位变量（保持与现有接口一致）
h_k0 = zeros(M,K); %#ok<NASGU>
f_k0 = zeros(N,K); %#ok<NASGU>
G0   = zeros(N,M); %#ok<NASGU>

%% ========== 轨迹记录 ==========
x_slot = zeros(Nslot+1,1);
y_slot = zeros(Nslot+1,1);
x_slot(1) = x_start;
y_slot(1) = y_start;
wsr_slot = zeros(Nslot,1);

%% ========== 局部 SCA 参数 ==========
inner_max_iter = 8;     % 每个时隙内 AO 次数
mu0 = 0.5;              % 保留接口（步长主控在 step_ratio）
grad_delta = 1.0;       % 数值梯度差分步长（米）
tol_pos = 1e-3;         % 位置收敛阈值
max_backtrack = 8;      % WSR 下降时最多减半次数
min_move_tol = 1e-4;    % 最小有效移动阈值

%% ========== 参数打包 ==========
params = struct();
params.K = K;
params.M = M;
params.N = N;
params.Ps_max = Ps_max;
params.Pr_max = Pr_max;
params.sigma2 = sigma2;
params.sigmar2 = sigmar2;
params.eta_k = eta_k;

params.bs_pos = bs_pos;
params.user_pos = user_pos;
params.h_const = h_RIS;

params.lambda = lambda;
params.gamma_reflect = gamma_reflect;
params.K_rician = K_rician;
params.scale_G = scale_G;
params.scale_f = scale_f;

params.v_max = v_max;
params.dt = dt;
params.x_min = x_min;
params.x_max = x_max;
params.y_min = y_min;
params.y_max = y_max;
params.step_ratio = 1.0;  % 归一化梯度步进系数（相对 v_max*dt）

%% =========================================
% 逐时隙联合优化
% ==========================================
for n = 1:Nslot
    fprintf('\n----------------------------------------\n');
    fprintf('Slot %d / %d\n', n, Nslot);

    q_curr = [x_slot(n), y_slot(n)];
    q_next = q_curr;

    fprintf('当前 USV 位置 = [%.2f, %.2f]\n', q_curr(1), q_curr(2));

    for it = 1:inner_max_iter
        % 1) 固定候选位置，构建单时隙信道
        ch = build_slot_channels_single(params, q_next);

        % 与主流程一致的直达链路衰减
        h_for_opt = 0.3 * ch.h_k;
        f_for_opt = ch.f_k;
        G_for_opt = ch.G;

        sigma2_eff  = sigma2  * (scale_G^2 * scale_f^2);
        sigmar2_eff = sigmar2 * (scale_G^2);
        Pr_max_eff  = Pr_max  * (scale_G^2);

        % 2) 固定位置，优化通信变量
        [W_new, Theta_new, Rsum_tmp, ~, ~] = active_RIS_precoding( ...
            M, K, N, Ps_max, Pr_max_eff, sigma2_eff, sigmar2_eff, eta_k, ...
            Theta_cur, W_cur, h_for_opt, f_for_opt, G_for_opt);

        % 3) 固定通信变量，数值计算位置梯度
        grad = numerical_gradient_slot_fixed_comm(q_next, W_new, Theta_new, params, grad_delta);

        % 4) 归一化梯度步进，得到候选点
        q_prop = trajectory_update_slot_localsca_from_grad(q_curr, q_next, grad, mu0, params);

        % 5) 单调性保护：若 WSR 下降，步长减半回溯
        wsr_ref = eval_slot_wsr_fixed_comm(q_next, W_new, Theta_new, params);
        dir_vec = q_prop - q_next;
        alpha = 1.0;
        q_new = q_prop;
        wsr_new = eval_slot_wsr_fixed_comm(q_new, W_new, Theta_new, params);

        bt = 0;
        while (wsr_new + 1e-12 < wsr_ref) && (bt < max_backtrack)
            alpha = alpha / 2;
            q_new = q_next + alpha * dir_vec;
            wsr_new = eval_slot_wsr_fixed_comm(q_new, W_new, Theta_new, params);
            bt = bt + 1;
        end

        if wsr_new + 1e-12 < wsr_ref
            q_new = q_next;
            wsr_new = wsr_ref;
            alpha = 0;
        end

        move_dist = norm(q_new - q_next);
        fprintf('  inner iter %d: Rsum=%.6f, grad=[%.4e, %.4e], move=%.4f m, alpha=%.3f, WSR: %.6f -> %.6f\n', ...
            it, Rsum_tmp, grad(1), grad(2), move_dist, alpha, wsr_ref, wsr_new);

        q_next = q_new;
        W_cur = W_new;
        Theta_cur = Theta_new;

        if move_dist < max(tol_pos, min_move_tol)
            break;
        end
    end

    % 更新到下一时隙位置
    x_slot(n+1) = q_next(1);
    y_slot(n+1) = q_next(2);

    % 用最终位置评估当前时隙 WSR
    wsr_slot(n) = eval_slot_wsr_fixed_comm(q_next, W_cur, Theta_cur, params);

    fprintf('Slot %d 完成: q_next = [%.2f, %.2f], WSR = %.6f\n', ...
        n, q_next(1), q_next(2), wsr_slot(n));
end

%% ========== 结果统计 ==========
[wsr_max, best_slot] = max(wsr_slot);
fprintf('\n========== 优化完成 ==========\n');
fprintf('平均 WSR = %.6f bps/Hz\n', mean(wsr_slot));
fprintf('最大 WSR = %.6f bps/Hz, 出现在 slot %d\n', wsr_max, best_slot);
fprintf('最终位置 = [%.2f, %.2f]\n', x_slot(end), y_slot(end));

%% ========== 绘图 ==========
figure('Position', [100, 100, 1200, 480]);

subplot(1,2,1); hold on; grid on; box on;
plot(bs_pos(1), bs_pos(2), 'ks', 'MarkerSize', 11, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
scatter(user_pos(:,1), user_pos(:,2), 55, 'b', 'filled');

th = linspace(0, 2*pi, 100);
plot(user_X_center + user_R*cos(th), user_Y_center + user_R*sin(th), '--', 'Color', [0.6 0.6 0.6]);

plot(x_slot, y_slot, '-r', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 10);
plot(x_start, y_start, 'mo', 'MarkerSize', 9, 'MarkerFaceColor', 'm');

step_mark = 5;
for t = 1:step_mark:(Nslot+1)
    text(x_slot(t), y_slot(t)+12, sprintf('%d', t-1), 'FontSize', 8, 'Color', [0.45 0 0.45]);
end

xlabel('x (m)');
ylabel('y (m)');
title('USV逐时隙轨迹（局部SCA）');
legend({'BS', 'Users', 'User Region', 'USV Trajectory', 'Initial Position'}, 'Location', 'best');
xlim([x_min-100, x_max+100]);
ylim([y_min-100, y_max+100]);

subplot(1,2,2); hold on; grid on; box on;
plot(1:Nslot, wsr_slot, '-r', 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
plot(best_slot, wsr_max, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
text(best_slot, wsr_max + 0.02*max(1,wsr_max), sprintf('max=%.3f', wsr_max), ...
    'HorizontalAlignment', 'center', 'FontSize', 10);

xlabel('Slot Index');
ylabel('WSR (bps/Hz)');
title('逐时隙 WSR');

save('result_1000m_slotwise_localsca.mat', ...
    'x_slot', 'y_slot', 'wsr_slot', 'params', 'user_pos', 'bs_pos', 'W_cur', 'Theta_cur');
