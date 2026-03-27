function Phi = Interference_cancellation_multiuser(theta,H_k,M,Upper_M,Iteration)

N = length(theta);

K = size(H_k,1);

Phi_opt=theta;

Phi = Phi_opt;

Phi_2 = Phi;

A = zeros(N,N);
b = zeros(N,1);

for a = 1:Iteration

    A = 0*A;
    b = 0*b;

    for k = 1:K
        H_k_temp = reshape(H_k(k,:,:),N,N);        
        temp = eye(N)+diag(Phi_2)*H_k_temp';
        A = A + (temp)'*(temp);
        b = b + (temp)'*Phi_opt;
    end
    
    Phi = ((A+M*eye(N))^(-1))*(M*Phi_2+b);

    A = 0*A;
    b = 0*b;

    for k = 1:K
        H_k_temp = reshape(H_k(k,:,:),N,N);
        temp = diag(H_k_temp'*Phi);
        A = A + (temp)'*(temp);
        b = b + (temp)'*(Phi_opt-Phi);
    end    

    Phi_2 = ((A+M*eye(N))^(-1))*(b+M*Phi);

    M = M*1.1;
    M = min(M,Upper_M);
     
end

end

