function PL = getPathLoss(dis, ht, hr, lambda)
% getPathLoss (modified version)
% 路径损耗采用你修改后的简化海上双射线公式：
%   PL(d) = (ht*hr)^2 / (4 * d^4)
% 说明：参数 lambda 在该版本中不参与计算，保留仅为接口兼容。

% 避免分母为 0
d = max(dis, eps);

% 路损简化公式
PL = (ht*hr)^2 / (4 * d^4);

end