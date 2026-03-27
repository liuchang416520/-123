function grad = numerical_gradient_slot_fixed_comm(q_xy, W, Theta, params, delta)
% 用中心差分计算当前位置处的 WSR 数值梯度
% grad = [dF/dx, dF/dy]

if nargin < 5
    delta = 1.0;
end

x = q_xy(1);
y = q_xy(2);

% x方向差分
F_x_plus  = eval_slot_wsr_fixed_comm([x + delta, y], W, Theta, params);
F_x_minus = eval_slot_wsr_fixed_comm([x - delta, y], W, Theta, params);
gx = (F_x_plus - F_x_minus) / (2 * delta);

% y方向差分
F_y_plus  = eval_slot_wsr_fixed_comm([x, y + delta], W, Theta, params);
F_y_minus = eval_slot_wsr_fixed_comm([x, y - delta], W, Theta, params);
gy = (F_y_plus - F_y_minus) / (2 * delta);

grad = [gx, gy];
end