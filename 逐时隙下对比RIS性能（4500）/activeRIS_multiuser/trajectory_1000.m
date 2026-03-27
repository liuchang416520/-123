function [x_new, y_new] = trajectory_update_sca(x_traj, y_traj, x_BR_0, x_RU_0, ...
    user_pos, bs_pos, dt, v_max, y_min, y_max, ...
    E1, E2, E3, tau_power, kappa_power, h_const, Ps_max, Pr_max, varargin)

    Nslot = numel(x_traj);
    K = size(user_pos,1);
    if ~isempty(varargin), lambda = varargin{1}; else, lambda = 0.06; end
    fix_endpoint = true;

    % 梯度权重（让算法更看重用户端的拉力）
    bs_weight = 0.3;
    ru_weight = 1.0;

    ht_BS = bs_pos(3); hr_RIS = h_const; hr_User = user_pos(1, 3);
    d_break_BR = (4 * ht_BS * hr_RIS) / lambda;
    d_break_RU = (4 * hr_RIS * hr_User) / lambda;

    U2 = zeros(K, Nslot);
    U3 = zeros(1, Nslot);

    for n = 1:Nslot
        x_br = max(x_BR_0(n), 1.0);
        if x_br <= d_break_BR, kappa_power_eff = 2; kappa_amp_eff = 1; else, kappa_power_eff = 4; kappa_amp_eff = 2; end
        
        for k = 1:K
            x_ru = max(x_RU_0(k,n), 1.0);
            if x_ru <= d_break_RU, tau_power_eff = 2; tau_amp_eff = 1; else, tau_power_eff = 4; tau_amp_eff = 2; end
            
            E1_nk = E1(k,n); E2_nk = E2(k,n); E3_nk = E3(k,n);
            term1 = - (tau_amp_eff/2) * (x_ru^(-tau_amp_eff/2 - 1)) * (x_br^(-kappa_amp_eff/2)) * E1_nk; 
            term2 =      tau_power_eff * (x_ru^(-tau_power_eff - 1)) * (x_br^(-kappa_power_eff)) * E2_nk;
            term3 =      tau_power_eff * (x_ru^(-tau_power_eff - 1)) * E3_nk; 
            U2(k,n) = real(term1 + term2 + term3);
        end

        term4_acc = 0; term5_acc = 0;
        for k = 1:K
            x_ru = max(x_RU_0(k,n), 1.0);
            E1_nk = E1(k,n); E2_nk = E2(k,n);
            if x_ru <= d_break_RU, tau_power_eff = 2; tau_amp_eff = 1; else, tau_power_eff = 4; tau_amp_eff = 2; end
            term4_acc = term4_acc + ( - (kappa_amp_eff/2) * (x_br^(-kappa_amp_eff/2 - 1)) * (x_ru^(-tau_amp_eff/2)) * E1_nk );
            term5_acc = term5_acc + (      kappa_power_eff * (x_br^(-kappa_power_eff - 1)) * (x_ru^(-tau_power_eff))   * E2_nk );
        end
        U3(n) = real(term4_acc + term5_acc);
    end

    U2(~isfinite(U2)) = 0; U3(~isfinite(U3)) = 0;

    grad_X = zeros(Nslot, 1);
    grad_Y = zeros(Nslot, 1);
    for n = 1:Nslot
        dx_br = x_traj(n) - bs_pos(1); dy_br = y_traj(n) - bs_pos(2);
        dist_br = max(x_BR_0(n), 1.0);
        gx_bs = bs_weight * U3(n) * (dx_br / dist_br);
        gy_bs = bs_weight * U3(n) * (dy_br / dist_br);
        
        gx_ru = 0; gy_ru = 0;
        for k = 1:K
            dx_ru = x_traj(n) - user_pos(k, 1); dy_ru = y_traj(n) - user_pos(k, 2);
            dist_ru = max(x_RU_0(k,n), 1.0);
            gx_ru = gx_ru + ru_weight * U2(k,n) * (dx_ru / dist_ru);
            gy_ru = gy_ru + ru_weight * U2(k,n) * (dy_ru / dist_ru);
        end
        grad_X(n) = gx_bs + gx_ru; grad_Y(n) = gy_bs + gy_ru;
    end
    
    grad_norm = sqrt(grad_X.^2 + grad_Y.^2);
    mean_grad = mean(grad_norm(grad_norm > 0));
    if mean_grad > 0
        % 【重点】：完全释放梯度，不再乘以0.4削弱
        grad_X = grad_X / mean_grad; grad_Y = grad_Y / mean_grad;
    end

    % 缩小场景后的信赖域步长 (约为单时隙最大位移 30 的量级)
    trust_step = 20;  

    try
        cvx_begin quiet
            cvx_solver sdpt3
            variables x_var(Nslot) y_var(Nslot)
            % 【重点】：最干净的 SCA 目标函数，只顺着梯度走，没有 alpha_ref 惩罚！
            obj = sum( grad_X .* (x_var - x_traj) + grad_Y .* (y_var - y_traj) );
            maximize(obj)
            subject to
                y_var >= y_min; y_var <= y_max;
                x_var(1) == x_traj(1); y_var(1) == y_traj(1);
                if fix_endpoint, x_var(Nslot) == x_traj(Nslot); y_var(Nslot) == y_traj(Nslot); end
                for n = 1:(Nslot-1), norm([x_var(n+1) - x_var(n); y_var(n+1) - y_var(n)]) <= v_max*dt; end
                for n = 1:Nslot, norm([x_var(n) - x_traj(n); y_var(n) - y_traj(n)]) <= trust_step; end
        cvx_end
        x_new = x_var; y_new = y_var;
    catch
        % 回退机制也只用真实梯度
        step_fb = 10;
        dir_x = grad_X; dir_y = grad_Y;
        gnorm = sqrt(dir_x.^2 + dir_y.^2) + 1e-12;
        x_new = x_traj + step_fb * (dir_x ./ gnorm);
        y_new = min(max(y_traj + step_fb * (dir_y ./ gnorm), y_min), y_max);
    end
end
