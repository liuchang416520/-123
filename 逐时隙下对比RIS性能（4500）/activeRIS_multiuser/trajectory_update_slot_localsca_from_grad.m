function q_new = trajectory_update_slot_localsca_from_grad(q_curr, q_ref, grad, mu, params)
% 单时隙局部 SCA 位置更新（按归一化梯度方向步进）

% 保留 mu 接口兼容性（当前步长主要由 step_ratio 控制）
if nargin < 4
    mu = 1.0;
end

r_max = params.v_max * params.dt;
g_norm = norm(grad);

% 默认步长：可达半径的 100%，可通过 params.step_ratio 调整
step_ratio = 1.0;
if isfield(params, 'step_ratio') && ~isempty(params.step_ratio)
    step_ratio = params.step_ratio;
end
step_ratio = max(step_ratio, 0);
step_len = min(step_ratio * r_max, r_max);

% 梯度过小时不移动
if g_norm < 1e-12 || step_len <= 0
    q_try = q_ref;
else
    q_try = q_ref + step_len * (grad / g_norm);
end

% 1) 投影到速度圆盘
vec = q_try - q_curr;
dist = norm(vec);
if dist > r_max
    vec = vec / dist * r_max;
    q_try = q_curr + vec;
end

% 2) 投影到边界
q_try(1) = min(max(q_try(1), params.x_min), params.x_max);
q_try(2) = min(max(q_try(2), params.y_min), params.y_max);

q_new = q_try;
end
