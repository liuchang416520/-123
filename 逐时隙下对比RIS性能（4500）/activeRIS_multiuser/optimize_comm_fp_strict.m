function [comm, E1, E2, E3, sum_rate_slot, aux] = optimize_comm_fp_strict(K, M, N_RIS, Ps_max, sigma2, sigmar2, eta_k, channels, is_active, Pr_max, varargin)
%OPTIMIZE_COMM_FP_STRICT
% Copy of optimize_comm_fp with two additions:
%   1) returns raw E1/E2/E3 before normalization/clipping;
%   2) returns aux struct for trajectory surrogate construction.

use_slot_warm_start = true;
Nslot = numel(channels);
comm   = repmat(struct('W', [], 'Theta', [], 'w_k', [], 'H_eff', [], ...
                       'P_bs', [], 'P_ris', [], 'Ps_max', [], 'Pr_max', []), Nslot, 1);
E1 = zeros(K, Nslot);
E2 = zeros(K, Nslot);
E3 = zeros(K, Nslot);
sum_rate_slot = zeros(Nslot, 1);

Last_Theta = diag(exp(1j * 2 * pi * rand(N_RIS, 1)));
Last_W     = exp(1j * 2 * pi * rand(K * M, 1)) * sqrt(Ps_max / max(K * M, 1));

traj_cfg = struct();
is_passive = false;

for n = 1:Nslot
    h_k = channels(n).h_k;
    f_k = channels(n).f_k;
    G   = channels(n).G;

    h_k_bar = channels(n).h_k_bar;
    f_k_bar = channels(n).f_k_bar;
    G_bar   = channels(n).G_bar;

    Theta_init = diag(exp(1j * 2 * pi * rand(N_RIS, 1)));
    W_init     = exp(1j * 2 * pi * rand(K * M, 1)) * sqrt(Ps_max / max(K * M, 1));

    ext_init_ok = false;
    init_slot = struct();
    traj_cfg = struct();

    if ~isempty(varargin) && numel(varargin) >= 1 && isstruct(varargin{1})
        init_slot = varargin{1};
        if isfield(init_slot, 'Theta') && isfield(init_slot, 'W') ...
                && ~isempty(init_slot.Theta) && ~isempty(init_slot.W)
            Theta_init = init_slot.Theta;
            W_init     = init_slot.W;
            ext_init_ok = true;
        end
    end

    if ~isempty(varargin) && numel(varargin) >= 2 && isstruct(varargin{2})
        traj_cfg = varargin{2};
    end

    is_passive = isfield(traj_cfg, 'ris_mode') && strcmpi(traj_cfg.ris_mode, 'passive');
    if is_passive
        ext_init_ok = false;
    end

    if ext_init_ok
        % keep external warm start
    elseif ~is_passive && isfield(traj_cfg, 'use_smart_init') && traj_cfg.use_smart_init
        if norm(h_k(1, :)) > 1e-12
            phase_anchor = angle(h_k(1, 1));
        else
            phase_anchor = 0;
        end

        theta_phase = phase_anchor - angle(f_k(1, :).') - angle(G(:, 1));
        Theta_tmp = diag(exp(1j * theta_phase));

        W0_mat = zeros(M, K);
        for k = 1:K
            h_eff = h_k(k, :) + f_k(k, :) * Theta_tmp * G;
            h_eff_norm = norm(h_eff);
            if h_eff_norm < 1e-12
                W0_mat(:, k) = ones(M, 1) / sqrt(M) * sqrt(Ps_max / K);
            else
                W0_mat(:, k) = (h_eff' / h_eff_norm) * sqrt(Ps_max / K);
            end
        end

        P_rx = 0;
        for k = 1:K
            P_rx = P_rx + norm(G * W0_mat(:, k))^2;
        end
        P_rx = P_rx + N_RIS * sigmar2;
        a_max = sqrt(Pr_max / max(P_rx, 1e-20));

        Theta_init = a_max * Theta_tmp;
        W_init     = W0_mat(:);
    elseif use_slot_warm_start && ~is_passive
        Theta_init = Last_Theta;
        W_init     = Last_W;
    end

    if nargin >= 10 && is_active
        [W, Theta, Rsum, Rho_k_n, eps_k_n] = active_RIS_precoding( ...
            M, K, N_RIS, Ps_max, Pr_max, sigma2, sigmar2, eta_k, ...
            Theta_init, W_init, h_k, f_k, G);
    elseif is_passive
        [W, Theta, Rsum, Rho_k_n, eps_k_n] = passive_RIS_precoding( ...
            M, K, N_RIS, Ps_max, sigma2, eta_k, ...
            Theta_init, W_init, h_k, f_k, G);
    else
        [W, Theta, Rsum] = random_RIS_precoding( ...
            M, K, N_RIS, Ps_max, sigma2, eta_k, ...
            Theta_init, W_init, h_k, f_k, G);
        Rho_k_n = ones(K, 1);
        eps_k_n = ones(K, 1);
    end

    Last_Theta = Theta;
    Last_W     = W;
    sum_rate_slot(n) = Rsum;
    w_k = w_k_generate(K, M, W);

    P_bs = real(W' * W);
    U = zeros(N_RIS, N_RIS);
    for kk = 1:K
        wtemp = reshape(w_k(kk, :), M, 1);
        U = U + diag(G * wtemp) * (diag(G * wtemp))';
    end
    U = U + sigmar2 * eye(N_RIS);
    theta_vec = Theta * ones(N_RIS, 1);
    P_ris = real(theta_vec' * U * theta_vec);

    H_k = H_k_generate(K, M, N_RIS, h_k, f_k, G, Theta);

    for k = 1:K
        f_bar_k_vec = reshape(f_k_bar(k, :), N_RIS, 1);
        cascade_bar_k_T = (f_bar_k_vec' * Theta * G_bar);
        h_bar_k_vec = reshape(h_k_bar(k, :), M, 1);

        chi_bar_k_n = 1 + Rho_k_n(k);
        rho_k_n = eps_k_n(k);
        wk = reshape(w_k(k, :), M, 1);

        term1_signal = real(conj(rho_k_n) * (cascade_bar_k_T * wk));
        term1 = 2 * sqrt(chi_bar_k_n) * term1_signal;

        interact_sum = 0;
        interf_sum_bar = 0;
        for j = 1:K
            wj = reshape(w_k(j, :), M, 1);
            direct_j = h_bar_k_vec' * wj;
            casc_j   = cascade_bar_k_T * wj;
            interact_sum = interact_sum + 2 * real(conj(direct_j) * casc_j);
            interf_sum_bar = interf_sum_bar + abs(casc_j)^2;
        end

        term2 = 2 * abs(rho_k_n)^2 * interact_sum;
        term3 = 2 * abs(rho_k_n)^2 * interf_sum_bar;

        E1(k, n) = term1 - term2;
        E2(k, n) = term3;
        E3(k, n) = abs(rho_k_n)^2 * norm(f_bar_k_vec' * Theta, 2)^2 * sigmar2;
    end

    comm(n).W = W;
    comm(n).Theta = Theta;
    comm(n).w_k = w_k;
    comm(n).H_eff = H_k;
    comm(n).P_bs = P_bs;
    comm(n).P_ris = P_ris;
    comm(n).Ps_max = Ps_max;
    comm(n).Pr_max = Pr_max;
end

% preserve raw coefficients before normalization/clipping
aux = struct();
aux.E1_raw = E1;
aux.E2_raw = E2;
aux.E3_raw = E3;
aux.Ps_max = Ps_max;
aux.Pr_max = Pr_max;
aux.sum_rate_slot = sum_rate_slot;
aux.comm = comm;

% normalized version for compatibility
E1n = E1; E2n = E2; E3n = E3;
e1_scale = prctile(abs(E1n(:)), 95);
e2_scale = prctile(abs(E2n(:)), 95);
e3_scale = prctile(abs(E3n(:)), 95);

if ~isfinite(e1_scale) || e1_scale < 1e-12, e1_scale = max(abs(E1n(:))); end
if ~isfinite(e2_scale) || e2_scale < 1e-12, e2_scale = max(abs(E2n(:))); end
if ~isfinite(e3_scale) || e3_scale < 1e-12, e3_scale = max(abs(E3n(:))); end
if ~isfinite(e1_scale) || e1_scale < 1e-12, e1_scale = 1.0; end
if ~isfinite(e2_scale) || e2_scale < 1e-12, e2_scale = 1.0; end
if ~isfinite(e3_scale) || e3_scale < 1e-12, e3_scale = 1.0; end

E1 = E1n / e1_scale;
E2 = E2n / e2_scale;
E3 = E3n / e3_scale;

e2_ratio_cap = 3.0;
if isfield(traj_cfg, 'e2_ratio_cap') && ~isempty(traj_cfg.e2_ratio_cap)
    e2_ratio_cap = traj_cfg.e2_ratio_cap;
end
med_e1 = median(abs(E1(:)));
med_e2 = median(abs(E2(:)));
e2_ratio = med_e2 / (med_e1 + 1e-12);
if isfinite(e2_ratio) && e2_ratio > e2_ratio_cap
    E2 = E2 * (e2_ratio_cap / e2_ratio);
end

e_clip = 6.0;
if isfield(traj_cfg, 'e_clip') && ~isempty(traj_cfg.e_clip)
    e_clip = traj_cfg.e_clip;
end
E1 = e_clip * tanh(E1 / e_clip);
E2 = e_clip * tanh(E2 / e_clip);
E3 = e_clip * tanh(E3 / e_clip);

aux.e1_scale = e1_scale;
aux.e2_scale = e2_scale;
aux.e3_scale = e3_scale;
aux.e2_ratio = e2_ratio;
end
