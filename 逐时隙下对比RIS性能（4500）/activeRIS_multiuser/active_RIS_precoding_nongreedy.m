function [W,Theta,Rsum,Rho_k,eps_k,varargout] = active_RIS_precoding_nongreedy( ...
    M,K,N,Ps_max,Pr_max,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G,varargin)
%ACTIVE_RIS_PRECODING_NONGREEDY
% Slot-by-slot but non-greedy: each slot solves communication and trajectory
% surrogate updates until the slot objective converges.
%
% Required traj_cfg additions:
%   q_final = [4500; -300]
%   fix_endpoint_last_slot = true
%   endpoint_reachability = true
%   eps_slot, MaxIter, MinIter
%
% Outputs are kept consistent with active_RIS_precoding in trajectory mode.

iteration = 30;
traj_mode = false;
traj_cfg = struct();

if ~isempty(varargin)
    traj_cfg = varargin{1};
    if isstruct(traj_cfg) && isfield(traj_cfg,'enable') && traj_cfg.enable
        traj_mode = true;
    end
end

log_level = 0;
if isfield(traj_cfg,'log_level') && ~isempty(traj_cfg.log_level)
    log_level = traj_cfg.log_level;
end

if ~traj_mode
    % fall back to original optimizer in fixed-position mode
    [W,Theta,Rsum,Rho_k,eps_k,varargout{1:nargout-5}] = active_RIS_precoding( ...
        M,K,N,Ps_max,Pr_max,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G,varargin{:});
    return;
end

if ~isfield(traj_cfg,'Nslot'), traj_cfg.Nslot = numel(traj_cfg.x_traj); end
if ~isfield(traj_cfg,'MaxIter'), traj_cfg.MaxIter = 12; end
if ~isfield(traj_cfg,'MinIter'), traj_cfg.MinIter = 3; end
if ~isfield(traj_cfg,'eps_slot'), traj_cfg.eps_slot = 1e-3; end
if ~isfield(traj_cfg,'y_min'), traj_cfg.y_min = -inf; end
if ~isfield(traj_cfg,'y_max'), traj_cfg.y_max = inf; end
if ~isfield(traj_cfg,'x_coast_min'), traj_cfg.x_coast_min = []; end
if ~isfield(traj_cfg,'no_direct_link'), traj_cfg.no_direct_link = false; end
if ~isfield(traj_cfg,'fix_endpoint_last_slot'), traj_cfg.fix_endpoint_last_slot = true; end
if ~isfield(traj_cfg,'endpoint_reachability'), traj_cfg.endpoint_reachability = true; end
if ~isfield(traj_cfg,'lambda_end'), traj_cfg.lambda_end = 1e-4; end
if ~isfield(traj_cfg,'q_final') || isempty(traj_cfg.q_final)
    traj_cfg.q_final = [traj_cfg.x_traj(end); traj_cfg.y_traj(end)];
end

scale_G = 1e3;
scale_f = 1e6;
if isfield(traj_cfg,'scale_G'), scale_G = traj_cfg.scale_G; end
if isfield(traj_cfg,'scale_f'), scale_f = traj_cfg.scale_f; end

Pr_max_eff   = Pr_max   * (scale_G^2);
sigmar2_eff  = sigmar2  * (scale_G^2);
sigma2_eff   = sigma2   * (scale_G^2 * scale_f^2);

x_opt = traj_cfg.x_traj(:);
y_opt = traj_cfg.y_traj(:);
sum_rate_slot = zeros(traj_cfg.Nslot, 1);
W_slots = cell(traj_cfg.Nslot, 1);
Theta_slots = cell(traj_cfg.Nslot, 1);
slot_info = cell(traj_cfg.Nslot, 1);

% Slot 1: communication only at prescribed start point.
[ch_init, ~, ~] = build_slot_channels(K, N, M, traj_cfg.bs_pos, traj_cfg.user_pos, ...
    x_opt(1), y_opt(1), traj_cfg.h_const, traj_cfg.lambda, traj_cfg.gamma_reflect, ...
    traj_cfg.K_rician, traj_cfg.no_direct_link, scale_G, scale_f);
ch_init(1).h_k = 0.3 * ch_init(1).h_k;
[comm_1, ~, ~, ~, sr_1, aux_1] = optimize_comm_fp_strict(K, M, N, Ps_max, sigma2_eff, sigmar2_eff, ...
    eta_k, ch_init(1), true, Pr_max_eff, struct(), traj_cfg);

sum_rate_slot(1) = sr_1;
W_slots{1} = comm_1.W;
Theta_slots{1} = comm_1.Theta;
slot_info{1} = aux_1;
prev_comm = comm_1;

if log_level >= 1
    fprintf('slot=1 fixed start: pos=(%.1f, %.1f), WSR=%.4f\n', x_opt(1), y_opt(1), sr_1);
end

for n = 2:traj_cfg.Nslot
    x_prev = x_opt(n-1);
    y_prev = y_opt(n-1);

    x_local = x_prev;
    y_local = y_prev;
    J_old = -inf;
    obj_hist = nan(traj_cfg.MaxIter,1);
    move_hist = nan(traj_cfg.MaxIter,1);

    for it = 1:traj_cfg.MaxIter
        [ch_local, x_BR_0, x_RU_0] = build_slot_channels(K, N, M, traj_cfg.bs_pos, traj_cfg.user_pos, ...
            x_local, y_local, traj_cfg.h_const, traj_cfg.lambda, traj_cfg.gamma_reflect, ...
            traj_cfg.K_rician, traj_cfg.no_direct_link, scale_G, scale_f);

        ch_local(1).h_k = 0.3 * ch_local(1).h_k;
        x_BR_0 = max(x_BR_0, 1e-3);
        x_RU_0 = max(x_RU_0, 1e-3);

        [comm_n, ~, ~, ~, sr_local, aux] = optimize_comm_fp_strict(K, M, N, Ps_max, sigma2_eff, sigmar2_eff, ...
            eta_k, ch_local(1), true, Pr_max_eff, prev_comm, traj_cfg);
        prev_comm = comm_n;

        traj_cfg.current_slot = n;
        [x_new, y_new, info] = trajectory_update_sca_strict( ...
            x_prev, y_prev, x_local, y_local, x_BR_0, x_RU_0, ...
            traj_cfg.user_pos, traj_cfg.bs_pos, traj_cfg.dt, traj_cfg.v_max, ...
            traj_cfg.y_min, traj_cfg.y_max, aux, traj_cfg);

        J_new = info.obj_lb;
        obj_hist(it) = J_new;
        move_hist(it) = norm([x_new - x_local; y_new - y_local]);
        rel_impr = inf;
        if isfinite(J_old)
            rel_impr = abs(J_new - J_old) / max(1, abs(J_old));
        end

        if log_level >= 1
            fprintf('  [slot=%d] it=%d, obj_lb=%.6e, rel=%.3e, move=%.3e -> (%.2f, %.2f)\n', ...
                n, it, J_new, rel_impr, move_hist(it), x_new, y_new);
        end

        x_local = x_new;
        y_local = y_new;

        if it >= traj_cfg.MinIter && rel_impr <= traj_cfg.eps_slot
            break;
        end
        J_old = J_new;
    end

    % one final communication update aligned with the converged slot location
    [ch_final, ~, ~] = build_slot_channels(K, N, M, traj_cfg.bs_pos, traj_cfg.user_pos, ...
        x_local, y_local, traj_cfg.h_const, traj_cfg.lambda, traj_cfg.gamma_reflect, ...
        traj_cfg.K_rician, traj_cfg.no_direct_link, scale_G, scale_f);
    ch_final(1).h_k = 0.3 * ch_final(1).h_k;
    [comm_n, ~, ~, ~, sr_final, aux_final] = optimize_comm_fp_strict(K, M, N, Ps_max, sigma2_eff, sigmar2_eff, ...
        eta_k, ch_final(1), true, Pr_max_eff, prev_comm, traj_cfg);
    prev_comm = comm_n;

    x_opt(n) = x_local;
    y_opt(n) = y_local;
    sum_rate_slot(n) = sr_final;
    W_slots{n} = comm_n.W;
    Theta_slots{n} = comm_n.Theta;
    slot_info{n} = struct('comm_aux', aux_final, 'obj_hist', obj_hist, 'move_hist', move_hist);

    if log_level >= 0
        fprintf('slot=%d converged: pos=(%.1f, %.1f), WSR=%.4f\n', n, x_opt(n), y_opt(n), sr_final);
    end
end

W = W_slots{end};
Theta = Theta_slots{end};
Rsum = sum_rate_slot(end);
Rho_k = [];
eps_k = [];

if nargout > 5
    varargout{1} = x_opt;
    varargout{2} = y_opt;
    varargout{3} = sum_rate_slot;
    if nargout > 8
        varargout{4} = slot_info;
    end
end
end
