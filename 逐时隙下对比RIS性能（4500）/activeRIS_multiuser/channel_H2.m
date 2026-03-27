function [H, H_bar] = channel_H2(N,M,dis,ht,hr,lambda,~,K)
% BS-User channel (修正版：使用断点路径损耗模型)
% 与 channel_F.m, channel_G.m 保持一致:
%   - 近场(d <= d_break): 功率衰减 d^-2 (幅度衰减 d^-1)
%   - 远场(d > d_break): 功率衰减 d^-4 (幅度衰减 d^-2)
%   - 断点距离: d_break = 4*ht*hr/lambda

d = max(dis, eps);
d_break = (4 * ht * hr) / lambda;

if d <= d_break
    % 近场: 幅度衰减 1/d
    const_amp = lambda / (4 * pi);
    dist_factor = 1 / d;
else
    % 远场: 幅度衰减 1/d^2
    const_amp = (ht * hr) / pi;
    dist_factor = 1 / d^2;
end

% Rician fading
H_LOS = exp(1j*2*pi*rand(N,M));
H_NLOS = (randn(N,M)+1j*randn(N,M))/sqrt(2);

% distance-independent component includes small-scale fading and constant amplitude
H_bar = const_amp * (sqrt(K/(K+1))*H_LOS + sqrt(1/(K+1))*H_NLOS);

% full channel with distance dependence
H = dist_factor * H_bar;

% Return row vector when M==1 to match callers expecting 1xN
if M == 1
    H = H.'; % 1 x N
    H_bar = H_bar.'; % 1 x N
end
end
