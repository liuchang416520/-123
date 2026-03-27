function [ h_k,f_k,G, h_k_bar,f_k_bar,G_bar] = Channel_generate2(K,N,M,Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser,ht_BS,h_RIS,hr_User,lambda,gamma_reflect,K_rician)
h_k=zeros(K,M);
f_k=zeros(K,N);
G=zeros(N,M);
h_k_bar=zeros(K,M);
f_k_bar=zeros(K,N);
G_bar=zeros(N,M);

for k=1:K
	[hk_full, hk_bar] = channel_H2(M,1,Dis_BStoUser(k),ht_BS,hr_User,lambda,gamma_reflect,K_rician);
	h_k(k,:) = hk_full;            
	h_k_bar(k,:) = hk_bar;
end


for k=1:K
	[f_full, f_bar] = channel_F(N,1,Dis_RIStoUser(k),h_RIS,hr_User,lambda,gamma_reflect,K_rician);
	f_k(k,:) = f_full;
	f_k_bar(k,:) = f_bar;
end

[G_full, Gbar] = channel_G(N,M,Dis_BStoRIS,ht_BS,h_RIS,lambda,gamma_reflect,K_rician);
G(:,:) = G_full;
G_bar(:,:) = Gbar;

end

