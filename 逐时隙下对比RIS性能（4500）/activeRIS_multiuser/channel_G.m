function [G, G_bar] = channel_G(N,M,dis,ht,hr,lambda,~,K,use_nlos)
% BS-RIS 反射前段信道 (使用统一的双射线模型)
%
% 关键修正：
%   - 基站高度 ht=40m，RIS高度 hr=3m
%   - 断点距离 d_break = 4*40*3/lambda ≈ 8.0 km (lambda≈0.06m)
%   - USV 正好在这个边缘活动 (6-10km)
%   - 分段模型能精确捕捉 USV 靠近基站时的能量红利
%   - 当 d < 8km 时使用 d^-2，当 d > 8km 时使用 d^-4

if nargin < 9 || isempty(use_nlos)
    use_nlos = true;
end

% 使用统一的路损模型
PL = getPathLoss_TwoRay(dis, ht, hr, lambda);

% Rician衰落（小尺度）
G_LOS = exp(1j*2*pi*rand(N,M));
if use_nlos
    G_NLOS = (randn(N,M)+1j*randn(N,M))/sqrt(2);
else
    G_NLOS = zeros(N,M);
end

% 分解为距离无关部分（用于SCA）
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

% 距离无关部分（小尺度衰落 + 常数幅度）
G_bar = const_amp * (sqrt(K/(K+1))*G_LOS + sqrt(1/(K+1))*G_NLOS);

% 完整信道 = 距离无关部分 * 距离因子
G = dist_factor * G_bar;
end
