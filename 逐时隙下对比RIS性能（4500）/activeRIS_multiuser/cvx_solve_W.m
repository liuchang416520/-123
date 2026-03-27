function W=cvx_solve_W(M,K,G,Theta,V,A,W,Ps_max,Pr_max,sigmar2)
% cvx_solve_W: 使用 CVX 求解基站波束赋形（FP 子问题最优解，QCQP）
%
% 优化问题:
%   minimize: real(W' * A * W) - 2 * real(V' * W)
%   subject to:
%     norm(W) <= sqrt(Ps_max)           (基站功率约束)
%     W' * D * W <= Pr_max - theta'*theta*sigmar2  (RIS功率约束)
%
% 输入:
%   M, K, G, Theta, V, A, W: 优化问题参数
%   Ps_max: 基站最大功率 (mW)
%   Pr_max: RIS最大功率 (mW)
%   sigmar2: RIS噪声功率
%
% 输出:
%   W: 优化后的波束赋形向量 (满足功率约束)

% 构建RIS功率约束矩阵
temp = (G'*Theta)*(Theta'*G);
theta = Theta*ones(size(Theta,1),1);
D = kron(eye(K), temp);
D = 0.5*(D+D');  % 确保D是Hermitian矩阵

% 确保A是Hermitian矩阵，并添加微小扰动保证正定性
A = 0.5*(A+A') + 1e-12*eye(size(A,1));

% RIS功率约束的右端项
P_ris_max = Pr_max - real(theta' * theta) * sigmar2;

% 使用CVX求解
try
    cvx_begin quiet
        cvx_solver sdpt3
        cvx_precision low   % 与 cvx_solve_theta 一致；'default' 在某些 CVX 版本可能触发 variable/char 异常
        
        variable W_var(M*K,1) complex
        
        % 目标函数: minimize real(W' * A * W) - 2 * real(V' * W)
        minimize(real(W_var' * A * W_var) - 2 * real(V' * W_var));
        
        subject to
            % 基站功率约束: norm(W) <= sqrt(Ps_max)
            norm(W_var) <= sqrt(Ps_max);
            
            % 主动RIS功率约束: W' * D * W <= Pr_max - theta'*theta*sigmar2
            real(W_var' * D * W_var) <= P_ris_max;
            
    cvx_end
    
    % 检查求解状态
    if strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Inaccurate/Solved')
        W = W_var;
        
        % 验证约束满足（仅用于调试）
        P_bs_check = real(W' * W);
        P_ris_check = real(W' * D * W);
        if P_bs_check > Ps_max * 1.01 || P_ris_check > P_ris_max * 1.01
            warning('cvx_solve_W: CVX求解后约束仍不满足！P_bs=%.2e (max=%.2e), P_ris=%.2e (max=%.2e)', ...
                P_bs_check, Ps_max, P_ris_check, P_ris_max);
        end
    else
        % CVX求解失败，使用输入值作为回退
        warning('cvx_solve_W: CVX求解失败 (状态=%s)，使用输入值', cvx_status);
        % 对输入W进行归一化，确保满足约束
        P_bs_input = real(W' * W);
        if P_bs_input > Ps_max
            W = W * sqrt(Ps_max / P_bs_input);
        end
    end
    
catch ME
    % CVX异常，使用输入值作为回退
    warning('cvx_solve_W: CVX异常 (%s)，使用输入值', ME.message);
    % 对输入W进行归一化，确保满足约束
    P_bs_input = real(W' * W);
    if P_bs_input > Ps_max
        W = W * sqrt(Ps_max / P_bs_input);
    end
end

end
