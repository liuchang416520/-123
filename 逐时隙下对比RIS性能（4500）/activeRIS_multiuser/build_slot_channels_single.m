function ch = build_slot_channels_single(params, q_xy)
% 单时隙信道生成（与 build_slot_channels 的缩放口径保持一致）

K = params.K;
N_RIS = params.N;
M = params.M;

bs_pos = params.bs_pos;
user_pos = params.user_pos;
h_const = params.h_const;
lambda = params.lambda;
gamma_reflect = params.gamma_reflect;
K_rician = params.K_rician;

scale_G = 1e3;
scale_f = 1e6;
if isfield(params, 'scale_G') && ~isempty(params.scale_G)
    scale_G = params.scale_G;
end
if isfield(params, 'scale_f') && ~isempty(params.scale_f)
    scale_f = params.scale_f;
end

no_direct_link = false;
if isfield(params, 'no_direct_link')
    no_direct_link = params.no_direct_link;
end

direct_link_attenuation = 0.354;

ris_pos = [q_xy(1), q_xy(2), h_const];

d_BR = max(norm(bs_pos - ris_pos), 1.0);

d_RU = zeros(K,1);
d_BU = zeros(K,1);
for k = 1:K
    d_RU(k) = max(norm(ris_pos - user_pos(k,:)), 1.0);
    d_BU(k) = max(norm(bs_pos - user_pos(k,:)), 1.0);
end

[h_k, f_k, G, h_k_bar, f_k_bar, G_bar] = Channel_generate2( ...
    K, N_RIS, M, d_BR, d_BU, d_RU, ...
    bs_pos(3), h_const, user_pos(1,3), ...
    lambda, gamma_reflect, K_rician);

% 直达链路控制
if no_direct_link
    h_k = zeros(size(h_k));
    h_k_bar = zeros(size(h_k_bar));
else
    h_k = h_k * direct_link_attenuation;
    h_k_bar = h_k_bar * direct_link_attenuation;
end

% 与主流程一致的分离缩放
G_scaled = G * scale_G;
f_k_scaled = f_k * scale_f;
h_k_scaled = h_k * (scale_G * scale_f);
G_bar_scaled = G_bar * scale_G;
f_k_bar_scaled = f_k_bar * scale_f;
h_k_bar_scaled = h_k_bar * (scale_G * scale_f);

ch = struct();
ch.h_k = h_k_scaled;
ch.f_k = f_k_scaled;
ch.G = G_scaled;
ch.h_k_bar = h_k_bar_scaled;
ch.f_k_bar = f_k_bar_scaled;
ch.G_bar = G_bar_scaled;
ch.d_BR = d_BR;
ch.d_RU = d_RU;
ch.d_BU = d_BU;
end
