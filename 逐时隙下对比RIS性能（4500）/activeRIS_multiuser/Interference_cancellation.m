function [cost, Phi] = Interference_cancellation(theta,H,M,Upper_M,Iteration)

N = length(theta);

Phi_opt=theta;

Phi = Phi_opt;
Phi_2 = Phi;

% H = SI_factor*1/sqrt(2)*(rand(N,N)+1j*rand(N,N));

%cost = zeros(Iteration,1);   
%cost =zeros;

for a = 1:Iteration
    % cost = norm(Phi+diag(Phi_2)*H'*Phi-Phi_opt)^2+M*norm(Phi-Phi_2)^2;
    
    A = eye(N)+diag(Phi_2)*H';
    b = -Phi_opt;
    
    Phi = ((A'*A+M*eye(N))^(-1))*(M*Phi_2-A'*b);
    
    cost = norm(Phi+diag(Phi_2)*H'*Phi-Phi_opt)^2+M*norm(Phi-Phi_2)^2;
    
    A = diag(H'*Phi);
    b = Phi-Phi_opt;
    
    Phi_2 = ((A'*A+M*eye(N))^(-1))*(M*Phi-A'*b);  
    
    M = M*1.1;
    M = min(M,Upper_M);
    
    if a>10 && norm(Phi-Phi_2)^2<0.001 && norm(Phi+diag(Phi)*H'*Phi-theta)^2<0.01
        break;
    end

end


end

