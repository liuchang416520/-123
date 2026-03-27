function [x_new, y_new, info] = trajectory_update_sca_strict( ...
    x_prev, y_prev, x_local, y_local, x_BR_0, x_RU_0, ...
    user_pos, bs_pos, dt, v_max, y_min, y_max, aux, traj_cfg)
%TRAJECTORY_UPDATE_SCA_STRICT
% Non-greedy one-slot trajectory update with objective-based convergence.
% It solves a convex surrogate subproblem around the current local point.
%
% Required fields in aux:
%   E1_raw, E2_raw, E3_raw
%
% Recommended fields in traj_cfg:
%   lambda, h_const, x_coast_min, trust_step, L_sca,
%   q_final, lambda_end, safe_user_r,
%   endpoint_reachability, fix_endpoint_last_slot,
%   current_slot, Nslot
%
% Output:
%   info.obj_lb    surrogate objective value
%   info.grad      gradient at local point
%   info.f0        proxy objective at local point
%   info.status    cvx status

ql = [x_local; y_local];
qp = [x_prev;  y_prev];

lambda = get_cfg(traj_cfg, 'lambda', 0.06);
h_const = get_cfg(traj_cfg, 'h_const', bs_pos(3));
trust_step = get_cfg(traj_cfg, 'trust_step', 0.8 * v_max * dt);
L_sca = get_cfg(traj_cfg, 'L_sca', 1.0);
x_coast_min = get_cfg(traj_cfg, 'x_coast_min', -inf);

% hard terminal data
qF = [];
if isfield(traj_cfg, 'q_final') && ~isempty(traj_cfg.q_final)
    qF = traj_cfg.q_final(:);
end
current_slot = get_cfg(traj_cfg, 'current_slot', []);
Nslot = get_cfg(traj_cfg, 'Nslot', []);
fix_endpoint_last_slot = get_cfg(traj_cfg, 'fix_endpoint_last_slot', true);
endpoint_reachability = get_cfg(traj_cfg, 'endpoint_reachability', true);
lambda_end = get_cfg(traj_cfg, 'lambda_end', 0.0);

% optional user safety radii
safe_user_r = [];
if isfield(traj_cfg, 'safe_user_r') && ~isempty(traj_cfg.safe_user_r)
    safe_user_r = traj_cfg.safe_user_r(:);
end

% obstacle fields are optional
obs_centers = [];
obs_r = [];
if isfield(traj_cfg, 'obs_centers') && ~isempty(traj_cfg.obs_centers)
    obs_centers = traj_cfg.obs_centers;
end
if isfield(traj_cfg, 'obs_r') && ~isempty(traj_cfg.obs_r)
    obs_r = traj_cfg.obs_r(:);
end

% evaluate proxy and gradient
[f0, grad] = eval_proxy_and_grad(ql, x_BR_0, x_RU_0, user_pos, bs_pos, aux, lambda, h_const);

grad = grad(:);
info = struct('obj_lb', -inf, 'grad', grad, 'f0', f0, 'status', 'not_run');

% build terminal feasibility radius
remaining_step = [];
if ~isempty(current_slot) && ~isempty(Nslot)
    remaining_step = max(Nslot - current_slot, 0) * v_max * dt;
end

% solve surrogate subproblem
try
    cvx_begin quiet
        cvx_solver sdpt3
        variables x_var y_var
        q = [x_var; y_var];

        % proximal lower surrogate of proxy objective
        obj = grad' * (q - ql) - 0.5 * L_sca * sum_square(q - ql);

        % soft terminal attraction for non-final slots
        if ~isempty(qF) && lambda_end > 0
            obj = obj - lambda_end * sum_square(q - qF);
        end

        maximize(obj)
        subject to
            % boundaries
            y_var >= y_min;
            y_var <= y_max;
            if isfinite(x_coast_min)
                x_var >= x_coast_min;
            end

            % one-slot mobility from previous committed point
            norm(q - qp) <= v_max * dt;

            % local trust region
            norm(q - ql) <= trust_step;

            % hard terminal equality on the last slot
            if fix_endpoint_last_slot && ~isempty(qF) && ~isempty(current_slot) && ~isempty(Nslot) && current_slot == Nslot
                x_var == qF(1);
                y_var == qF(2);
            end

            % remaining-horizon terminal reachability
            if endpoint_reachability && ~isempty(qF) && ~isempty(remaining_step) && current_slot < Nslot
                norm(q - qF) <= remaining_step;
            end

            % user safety constraints via first-order inner approximation
            if ~isempty(safe_user_r)
                for k = 1:size(user_pos,1)
                    uk = user_pos(k,1:2).';
                    lhs_u = sum_square(ql - uk) + 2 * (ql - uk)' * (q - ql);
                    lhs_u >= safe_user_r(min(k, numel(safe_user_r)))^2;
                end
            end

            % obstacle avoidance via first-order inner approximation
            if ~isempty(obs_centers) && ~isempty(obs_r)
                for j = 1:size(obs_centers,1)
                    oj = obs_centers(j,1:2).';
                    lhs_o = sum_square(ql - oj) + 2 * (ql - oj)' * (q - ql);
                    lhs_o >= obs_r(min(j, numel(obs_r)))^2;
                end
            end
    cvx_end

    x_new = x_var;
    y_new = y_var;
    info.obj_lb = cvx_optval;
    info.status = cvx_status;

    if ~isfinite(x_new) || ~isfinite(y_new)
        x_new = x_local;
        y_new = y_local;
        info.status = 'fallback_invalid';
    end
catch ME
    % fail-safe: keep current point, which preserves feasibility of slot-by-slot scheme
    x_new = x_local;
    y_new = y_local;
    info.status = ['fallback_exception: ' ME.identifier];
end

end

function [f0, grad] = eval_proxy_and_grad(q, x_BR_0, x_RU_0, user_pos, bs_pos, aux, lambda, h_const)
% Proxy objective:
%   sum_k [ E1_raw * phi_k(q) - E2_raw * psi_k(q) - gamma_damp * E3_raw * nu_k(q) ]
% where phi_k = d_BR^{-kappa_a/2} d_RU^{-tau_a/2},
%       psi_k = d_BR^{-kappa_p}   d_RU^{-tau_p},
%       nu_k  = d_RU^{-tau_p}.

K = size(user_pos,1);
ht_BS   = bs_pos(3);
hr_RIS  = h_const;
hr_User = user_pos(1,3);

% breakpoint distances
if isempty(lambda) || lambda <= 0
    lambda = 0.06;
end
d_break_BR = (4 * ht_BS * hr_RIS) / lambda;
d_break_RU = (4 * hr_RIS * hr_User) / lambda;

x_br = max(x_BR_0, 1.0);
if x_br <= d_break_BR
    kappa_power_eff = 2;
    kappa_amp_eff   = 1;
else
    kappa_power_eff = 4;
    kappa_amp_eff   = 2;
end

E1 = aux.E1_raw(:);
E2 = aux.E2_raw(:);
E3 = aux.E3_raw(:);

Ps_max = get_struct(aux, 'Ps_max', 1.0);
Pr_max = get_struct(aux, 'Pr_max', 0.0);
gamma_damping = min(10 * Pr_max / (Ps_max + 1e-6), 1.0);

f0 = 0;
grad = zeros(2,1);

% derivatives of d_BR wrt q
br_vec = q - bs_pos(1:2).';
if x_br < 1e-9
    dbr = [0;0];
else
    dbr = br_vec / x_br;
end

for k = 1:K
    x_ru = max(x_RU_0(k), 1.0);
    if x_ru <= d_break_RU
        tau_power_eff = 2;
        tau_amp_eff   = 1;
    else
        tau_power_eff = 4;
        tau_amp_eff   = 2;
    end

    uk = user_pos(k,1:2).';
    ru_vec = q - uk;
    dru = ru_vec / x_ru;

    phi = x_br^(-kappa_amp_eff/2) * x_ru^(-tau_amp_eff/2);
    psi = x_br^(-kappa_power_eff)  * x_ru^(-tau_power_eff);
    nu  = x_ru^(-tau_power_eff);

    f0 = f0 + E1(k) * phi - E2(k) * psi - gamma_damping * E3(k) * nu;

    % gradients
    dphi = phi * (-(kappa_amp_eff/2) * dbr / x_br - (tau_amp_eff/2) * dru / x_ru);
    dpsi = psi * (-(kappa_power_eff)   * dbr / x_br - (tau_power_eff) * dru / x_ru);
    dnu  = nu  * (-(tau_power_eff) * dru / x_ru);

    grad = grad + E1(k) * dphi - E2(k) * dpsi - gamma_damping * E3(k) * dnu;
end

end

function v = get_cfg(s, name, default_v)
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default_v;
end
end

function v = get_struct(s, name, default_v)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default_v;
end
end
