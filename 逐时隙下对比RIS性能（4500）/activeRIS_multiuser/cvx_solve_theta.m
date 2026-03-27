function theta = cvx_solve_theta(N,K,M,theta_old,nu,Lam,w_k,G,Pr_max,sigmar2)
% cvx_solve_theta
% 稳健版 theta 更新：不用 CVX，改为“闭式线性求解 + 投影 + 回退保护”
%
% 目标近似为：
%   maximize  -theta' * Lam * theta + 2*real(theta' * nu)
%             - beta_prox * ||theta - theta_old||^2
%   subject to theta' * U * theta <= Pr_max
%
% 处理方法：
% 1) 先构造 U，并做 Hermitian 化与正则；
% 2) 用闭式线性系统求无约束候选解；
% 3) 若候选解不满足功率约束，则按 U-范数投影缩放；
% 4) 若数值异常，则回退 theta_old。
%
% 这样可以避免 CVX/SDPT3 的 infeasible、RCOND=NaN、scaling error 等问题。

% ---------- 基本保护 ----------
theta_old = theta_old(:);
nu = nu(:);

if numel(theta_old) ~= N
    error('cvx_solve_theta: theta_old dimension mismatch.');
end
if numel(nu) ~= N
    error('cvx_solve_theta: nu dimension mismatch.');
end

% ---------- 构造 U ----------
U = zeros(N,N);
for k = 1:K
    w_k_temp = reshape(w_k(k,:), M, 1);
    gw = G * w_k_temp;              % N x 1
    Dk = diag(gw);
    U = U + Dk * Dk';
end
U = U + sigmar2 * eye(N);

% ---------- Hermitian 化 ----------
Lam = 0.5 * (Lam + Lam');
U   = 0.5 * (U   + U');

% ---------- 数值尺度 ----------
lam_scale = max(norm(Lam, 'fro'), 1e-12);
u_scale   = max(norm(U,   'fro'), 1e-12);
nu_scale  = max(norm(nu),         1e-12);
th_scale  = max(norm(theta_old),  1e-12);

% ---------- 正则化 ----------
% Lam 若太弱或半奇异，直接增强对角正则，避免线性系统病态
reg_lam = max(1e-8, 1e-4 * lam_scale);
reg_u   = max(1e-8, 1e-8 * u_scale);

Lam_reg = Lam + reg_lam * eye(N);
U_reg   = U   + reg_u   * eye(N);

% ---------- 近端项 ----------
% 近端项要足够强，防止 theta 一步骤降或乱跳
beta_prox = max(1e-4, 5e-2 * max(1, lam_scale));

% ---------- 右端项 ----------
rhs = nu + beta_prox * theta_old;

% ---------- 构造线性系统 ----------
% 无约束最优近似：
%   (Lam + beta I) theta = nu + beta theta_old
A = Lam_reg + beta_prox * eye(N);

% 再做一次 Hermitian 化，避免微小非对称引起数值问题
A = 0.5 * (A + A');

% ---------- 条件数保护 ----------
% 若 A 仍接近奇异，再增大对角加载
rcondA = rcond(A);
if ~isfinite(rcondA) || rcondA < 1e-12
    add_diag = max(1e-6, 1e-2 * trace(abs(A)) / max(N,1));
    A = A + add_diag * eye(N);
end

% ---------- 闭式求解 ----------
theta = theta_old;
solve_ok = true;

try
    theta_cand = A \ rhs;
catch
    solve_ok = false;
    theta_cand = theta_old;
end

% ---------- 数值有效性检查 ----------
if ~solve_ok || any(~isfinite(theta_cand))
    theta_cand = theta_old;
end

% ---------- 崩塌保护 ----------
old_norm = max(norm(theta_old), 1e-12);
cand_norm = norm(theta_cand);

% 候选解若太小，说明掉进零陷阱，直接回退
if cand_norm < 0.2 * old_norm
    theta_cand = theta_old;
end

% 候选解若暴涨太多，也做限制
if cand_norm > 5.0 * old_norm
    theta_cand = theta_cand * (5.0 * old_norm / max(cand_norm,1e-12));
end

% ---------- 投影到 RIS 功率约束 ----------
P_cand = real(theta_cand' * U_reg * theta_cand);

if ~isfinite(P_cand) || P_cand < 0
    theta_cand = theta_old;
    P_cand = real(theta_cand' * U_reg * theta_cand);
end

if P_cand > Pr_max && P_cand > 1e-20
    % 用 U-范数缩放到边界
    theta_cand = theta_cand * sqrt(Pr_max / P_cand);
end

% ---------- 最终平滑更新 ----------
% 再做一层保守混合，增强稳定性
alpha_local = 0.3;
theta = (1 - alpha_local) * theta_old + alpha_local * theta_cand;

% ---------- 最终功率投影 ----------
P_final = real(theta' * U_reg * theta);
if P_final > Pr_max && P_final > 1e-20
    theta = theta * sqrt(Pr_max / P_final);
end

% ---------- 最终兜底 ----------
if any(~isfinite(theta))
    theta = theta_old;
end

% 若旧解本身也异常，再随机相位兜底
if any(~isfinite(theta_old))
    theta = exp(1j * 2*pi * rand(N,1));
    P_rand = real(theta' * U_reg * theta);
    if P_rand > Pr_max && P_rand > 1e-20
        theta = theta * sqrt(Pr_max / P_rand);
    end
end

end