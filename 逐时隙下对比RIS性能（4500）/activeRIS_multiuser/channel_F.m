function [F, F_bar] = channel_F(N,M,dis,ht,hr,lambda,~,K)


% 使用统一的路损模型
PL = getPathLoss_TwoRay(dis, ht, hr, lambda);

% Rician衰落（小尺度）
F_LOS = exp(1j*2*pi*rand(N,M));
F_NLOS = (randn(N,M)+1j*randn(N,M))/sqrt(2);

% 分解为距离无关部分（用于SCA）
% sqrt(PL) = loss_amp，需要分解为 const_amp * dist_factor
% 对于近场: sqrt(PL) = lambda/(4*pi*d) = [lambda/(4*pi)] * [1/d]
% 对于远场: sqrt(PL) = h1*h2/(pi*d^2) = [h1*h2/pi] * [1/d^2]

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
F_bar = const_amp * (sqrt(K/(K+1))*F_LOS + sqrt(1/(K+1))*F_NLOS);

% 完整信道 = 距离无关部分 * 距离因子
F = dist_factor * F_bar;

% 返回行向量（M==1时）
if M == 1
    F = F.';
    F_bar = F_bar.';
end
end
