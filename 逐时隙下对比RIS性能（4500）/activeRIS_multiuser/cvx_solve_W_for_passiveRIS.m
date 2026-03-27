function W=cvx_solve_W_for_passiveRIS(M,K,G,Theta,V,A,W,Ps_max)
% 修复矩阵病态：归一化 + 正则化，避免 RCOND 接近奇异警告

% 确保 Hermitian
A = 0.5*(A+A');

% 归一化：将 A、V 缩放到合理量级，改善条件数
scale = norm(A,'fro');
if scale < 1e-20
    scale = 1;
end
A = A / scale;
V = V / scale;

% 正则化：添加小量单位阵，保证正定（过大会削弱信道自适应，导致WSR趋平）
delta = 1e-6;
A_reg = A + delta * eye(size(A,1));

% 使用 mldivide 替代 inv()，数值更稳定
W = A_reg \ V;

if real(W'*W) <= Ps_max
    return;
end

% 二分法求满足功率约束的 lambda
lambda_left = 0;
lambda_right = 10;

W_right = (A_reg + lambda_right*eye(size(A,1))) \ V;
while real(W_right'*W_right) > Ps_max
    lambda_right = 2 * lambda_right;
    W_right = (A_reg + lambda_right*eye(size(A,1))) \ V;
end

lambda_center = (lambda_left + lambda_right) / 2;
W_center = (A_reg + lambda_center*eye(size(A,1))) \ V;

while (lambda_right - lambda_left) / max(lambda_right, 1e-10) > 1e-7
    if real(W_center'*W_center) < Ps_max
        lambda_right = lambda_center;
    else
        lambda_left = lambda_center;
    end
    lambda_center = (lambda_left + lambda_right) / 2;
    W_center = (A_reg + lambda_center*eye(size(A,1))) \ V;
end

W = W_center;

end