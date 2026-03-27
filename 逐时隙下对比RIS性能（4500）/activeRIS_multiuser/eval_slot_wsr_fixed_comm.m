function Rsum = eval_slot_wsr_fixed_comm(q_xy, W, Theta, params)

% 在固定通信变量 (W, Theta) 下，评估位置 q_xy 对应的加权和速率

ch = build_slot_channels_single(params, q_xy);

K = params.K;
M = params.M;
eta_k = params.eta_k;

scale_G = params.scale_G;
scale_f = params.scale_f;

sigma2_eff  = params.sigma2  * (scale_G^2 * scale_f^2);
sigmar2_eff = params.sigmar2 * (scale_G^2);

h_k = 0.3 * ch.h_k;   % 与主流程一致的直达链路衰减
f_k = ch.f_k;
G   = ch.G;

w_k = split_beamformers(W, K, M);   % K x M

Rsum = 0;
for k = 1:K
    % Channel_generate2 输出维度：h_k(k,:), f_k(k,:)
    h_eq = h_k(k,:) + f_k(k,:) * Theta * G;    % 1 x M

    signal = abs(h_eq * w_k(k,:).')^2;

    interf = 0;
    for j = 1:K
        if j ~= k
            interf = interf + abs(h_eq * w_k(j,:).')^2;
        end
    end

    noise_ris = norm(f_k(k,:) * Theta, 2)^2 * sigmar2_eff;
    SINR_k = signal / (interf + noise_ris + sigma2_eff);

    Rsum = Rsum + eta_k(k) * log2(1 + max(real(SINR_k), 0));
end
end
