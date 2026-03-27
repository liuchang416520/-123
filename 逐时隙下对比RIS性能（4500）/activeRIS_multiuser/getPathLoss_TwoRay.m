function PL = getPathLoss_TwoRay(d, h1, h2, lambda)
% getPathLoss_TwoRay: 统一的双射线路径损耗模型
% 输入:
%   d      - 距离 (m)
%   h1     - 发射端高度 (m)
%   h2     - 接收端高度 (m)
%   lambda - 波长 (m)
% 输出:
%   PL     - 路径损耗增益 (线性，非dB)
%
% 物理原理：
%   - 近场(d <= d_break): 视距主导，使用自由空间模型 d^-2
%   - 远场(d > d_break): 干涉主导，使用双射线渐近模型 d^-4
%   - 断点: d_break = 4*h1*h2/lambda (第一菲涅尔区理论)

% 避免距离为0
d = max(d, eps);

% 计算链路特有的断点距离
d_break = (4 * h1 * h2) / lambda;

if d <= d_break
    % 近场/视距主导区: 自由空间传播模型 (功率衰减 d^-2)
    % PL = (lambda / (4*pi*d))^2
    PL = (lambda / (4 * pi * d))^2;
else
    % 远场/干涉主导区: 双射线渐近模型 (功率衰减 d^-4)
    % 在断点处确保连续性
    % 双射线公式: PL = (h1*h2)^2 / d^4
    % 为了连续，需要乘以一个修正系数
    % 在 d = d_break 处，两个模型应该相等:
    % (lambda/(4*pi*d_break))^2 = C * (h1*h2)^2 / d_break^4
    % 解出: C = (lambda/(4*pi))^2 * d_break^2 / (h1*h2)^2
    %        = (lambda/(4*pi))^2 * (4*h1*h2/lambda)^2 / (h1*h2)^2
    %        = (lambda/(4*pi))^2 * 16*(h1*h2)^2/lambda^2 / (h1*h2)^2
    %        = 16/(4*pi)^2 = 1/(pi^2)
    
    % 简化后的双射线模型（已保证连续）
    PL = (h1 * h2)^2 / (pi^2 * d^4);
end

end
