function [x_new, y_new] = trajectory_update_sca(x_prev, y_prev, x_local, y_local, x_BR_0, x_RU_0, ...
    user_pos, bs_pos, dt, v_max, y_min, y_max, ...
    E1, E2, E3, tau_power, kappa_power, h_const, Ps_max, Pr_max, varargin)

% 单时隙轨迹 SCA 更新器
% 改动重点：
% 1) 不再把梯度单位化，只做截断，保留强弱信息
% 2) 增加用户中心吸引项
% 3) 给 y 下边界加软 margin，减少贴边
% 4) 保留平台探索，但减小探索步长

K = size(user_pos,1);

if ~isempty(varargin)
    lambda = varargin{1};
else
    lambda = 0.06;
end

x_coast_min = -inf;
if numel(varargin) >= 2 && ~isempty(varargin{2})
    x_coast_min = varargin{2};
end

% =========================
% 断点计算
% =========================
ht_BS   = bs_pos(3);
hr_RIS  = h_const;
hr_User = user_pos(1,3);

d_break_BR = (4 * ht_BS * hr_RIS) / lambda;
d_break_RU = (4 * hr_RIS * hr_User) / lambda;

% =========================
% 有源 RIS 噪声阻尼
% =========================
if isempty(Ps_max) || isempty(Pr_max)
    gamma_damping = 1.0;
else
    gamma_damping = Pr_max / (Ps_max + 1e-6);
    gamma_damping = min(gamma_damping * 10, 1.0);
end

% =========================
% E2 软缩放，防止干扰项过强主导梯度
% =========================
e1_med = median(abs(E1(:)));
e2_med = median(abs(E2(:)));
e2_soft_cap = 3.0;
e2_soft_scale = 1.0;

if isfinite(e1_med) && isfinite(e2_med) && e1_med > 0 && e2_med > e2_soft_cap * e1_med
    e2_soft_scale = (e2_soft_cap * e1_med) / (e2_med + 1e-12);
end
E2_eff = E2 * e2_soft_scale;

% =========================
% 距离项
% =========================
U2 = zeros(K,1);

x_br = max(x_BR_0, 1.0);
if x_br <= d_break_BR
    kappa_power_eff = 2;
    kappa_amp_eff   = 1;
else
    kappa_power_eff = 4;
    kappa_amp_eff   = 2;
end

term4_acc = 0;
term5_acc = 0;

for k = 1:K
    x_ru = max(x_RU_0(k), 1.0);

    if x_ru <= d_break_RU
        tau_power_eff = 2;
        tau_amp_eff   = 1;
    else
        tau_power_eff = 4;
        tau_amp_eff   = 2;
    end

    E1_k = E1(k);
    E2_k = E2_eff(k);
    E3_k = E3(k);

    % 信号引力 / 干扰斥力 / 噪声斥力
    term1 = - (tau_amp_eff / 2) * (x_ru^(-tau_amp_eff/2 - 1)) * (x_br^(-kappa_amp_eff/2)) * E1_k;
    term2 =   tau_power_eff      * (x_ru^(-tau_power_eff   - 1)) * (x_br^(-kappa_power_eff)) * E2_k;
    term3 = gamma_damping * tau_power_eff * (x_ru^(-tau_power_eff - 1)) * E3_k;

    U2(k) = real(term1 + term2 + term3);

    term4_acc = term4_acc + (-(kappa_amp_eff / 2) * (x_br^(-kappa_amp_eff/2 - 1)) * (x_ru^(-tau_amp_eff/2)) * E1_k);
    term5_acc = term5_acc + (  kappa_power_eff     * (x_br^(-kappa_power_eff - 1)) * (x_ru^(-tau_power_eff))   * E2_k);
end

U3 = real(term4_acc + term5_acc);

U2(~isfinite(U2)) = 0;
U3(~isfinite(U3)) = 0;

% =========================
% 梯度投影
% =========================
dx_br = x_local - bs_pos(1);
dy_br = y_local - bs_pos(2);

gx_bs = U3 * (dx_br / x_br);
gy_bs = U3 * (dy_br / x_br);

gx_ru = 0;
gy_ru = 0;

for k = 1:K
    dx_ru = x_local - user_pos(k,1);
    dy_ru = y_local - user_pos(k,2);
    dist_ru = max(x_RU_0(k), 1.0);

    gx_ru = gx_ru + U2(k) * (dx_ru / dist_ru);
    gy_ru = gy_ru + U2(k) * (dy_ru / dist_ru);
end

grad_X = real(gx_bs + gx_ru);
grad_Y = real(gy_bs + gy_ru);

% =========================
% 用户中心方向
% =========================
user_center = mean(user_pos(:,1:2), 1);
dir_user = [user_center(1) - x_local; user_center(2) - y_local];
dir_user_norm = norm(dir_user);

if dir_user_norm > 1e-12
    dir_user = dir_user / dir_user_norm;
else
    dir_user = [1; 0];
end

% =========================
% 不再单位化梯度，只做截断
% =========================
g = [grad_X; grad_Y];
g_norm = norm(g);

if g_norm < 1e-12
    % 纯平台时，退化为用户中心方向
    grad_X = dir_user(1);
    grad_Y = dir_user(2);
else
    % 梯度截断而不是单位化
    g_cap = 1e3;
    scale_g = max(1.0, g_norm / g_cap);
    grad_X = grad_X / scale_g;
    grad_Y = grad_Y / scale_g;
end

% =========================
% trust region
% =========================
trust_step = 0.8 * (v_max * dt);

% =========================
% 目标中的用户吸引项
% =========================
mu_user = 5e-4;      % 用户中心吸引强度
mu_x_forward = 1e-6; % 向右推进偏置（保留但减弱）

% y 软 margin：避免长期贴 y_min
y_margin = 2.0;

try
    cvx_begin quiet
        cvx_solver sdpt3
        variables x_var y_var

        obj = grad_X * (x_var - x_local) ...
            + grad_Y * (y_var - y_local) ...
            + mu_user * (dir_user(1) * (x_var - x_local) + dir_user(2) * (y_var - y_local)) ...
            + mu_x_forward * max(0, user_center(1) - x_local) * 1e-3 * (x_var - x_prev);

        maximize(obj)

        subject to
            if isfinite(x_coast_min)
                x_var >= x_coast_min;
            end

            % y 边界 + 软 margin
            if y_local <= y_min + y_margin
                y_var >= y_min + y_margin;
            else
                y_var >= y_min;
            end
            y_var <= y_max;

            % 相对上一时隙的机动约束
            norm([x_var - x_prev; y_var - y_prev]) <= v_max * dt;

            % trust region
            norm([x_var - x_local; y_var - y_local]) <= trust_step;

            % 若用户区在右侧，则不允许 x 回退
            if mean(user_pos(:,1)) > x_local
                x_var >= x_local;
            end
    cvx_end

    x_new = x_var;
    y_new = y_var;

    if isnan(x_new) || isnan(y_new) || ~isfinite(x_new) || ~isfinite(y_new)
        x_new = x_local;
        y_new = y_local;
    end

    % =========================
    % anti-plateau exploration
    % =========================
    move_norm = norm([x_new - x_local; y_new - y_local]);
    if move_norm < 1e-3
        explore_step = 0.05 * v_max * dt;
        x_try = x_local + explore_step * dir_user(1);
        y_try = y_local + explore_step * dir_user(2);

        if y_local <= y_min + y_margin
            y_try = max(y_try, y_min + y_margin);
        end
        y_try = min(max(y_try, y_min), y_max);

        if isfinite(x_coast_min)
            x_try = max(x_try, x_coast_min);
        end

        if mean(user_pos(:,1)) > x_local
            x_try = max(x_try, x_local);
        end

        if norm([x_try - x_prev; y_try - y_prev]) <= v_max * dt + 1e-9
            x_new = x_try;
            y_new = y_try;
        end
    end

catch
    % CVX 失败时，给一个小的用户方向探索步，而不是完全原地不动
    explore_step = 0.03 * v_max * dt;

    x_new = x_local + explore_step * dir_user(1);
    y_new = y_local + explore_step * dir_user(2);

    if y_local <= y_min + y_margin
        y_new = max(y_new, y_min + y_margin);
    end
    y_new = min(max(y_new, y_min), y_max);

    if isfinite(x_coast_min)
        x_new = max(x_new, x_coast_min);
    end

    if mean(user_pos(:,1)) > x_local
        x_new = max(x_new, x_local);
    end

    if norm([x_new - x_prev; y_new - y_prev]) > v_max * dt + 1e-9
        x_new = x_local;
        y_new = y_local;
    end
end

end