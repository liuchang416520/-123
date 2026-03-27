function [H] = channel_H(N,M,dis,ht,hr,lambda,~,K)
% BS-User 直射信道 (使用统一的双射线模型)
% 
% 关键修正：
%   - 基站高度 ht=40m，用户高度 hr=8m
%   - 断点距离 d_break = 4*40*8/lambda ≈ 21.3 km (lambda≈0.06m)
%   - 在你的仿真范围内 (0-16km)，直射链路一直处于近场 d^-2 状态
%   - 之前错误地使用 d^-4，严重低估了直射信号强度
%   - 修正后，直射链路会变强，RIS的作用更多体现在"补盲"或"增强边缘用户"

% 使用统一的路损模型
PL = getPathLoss_TwoRay(dis, ht, hr, lambda);

% Rician衰落（小尺度）
H_LOS = exp(1j*2*pi*rand(N,M));
H_NLOS = (randn(N,M)+1j*randn(N,M))/sqrt(2);
H = sqrt(PL)*(sqrt(K/(K+1))*H_LOS + sqrt(1/(K+1))*H_NLOS);

% 返回行向量（M==1时）
if M == 1
    H = H.'; % 1 x N
end
end
