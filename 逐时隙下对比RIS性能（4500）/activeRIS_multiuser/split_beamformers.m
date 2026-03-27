function w_k = split_beamformers(W, K, M)
% 将 W 转换成 K x M，每一行对应一个用户的波束向量
%
% 兼容两种输入：
% 1) W 是 (K*M) x 1 向量（你原代码风格）
% 2) W 是 M x K 矩阵

if isvector(W)
    W = W(:);
    if length(W) ~= K*M
        error('W length mismatch: expected K*M.');
    end
    w_k = zeros(K, M);
    for k = 1:K
        idx = (k-1)*M + (1:M);
        w_k(k,:) = W(idx).';
    end
else
    [r,c] = size(W);
    if r == M && c == K
        w_k = W.';
    elseif r == K && c == M
        w_k = W;
    else
        error('Unsupported W size.');
    end
end
end