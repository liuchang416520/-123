tic

user_X=[100:50:600];
ExNumber=500; %2000

M=4;  %基站天线数
Ps_max=10;        %基站最大功率限制 mW
Pr_max=10;        %PA-RIS最大功率限制 mW

sigma2 = 10^(-10);

sigmar2=sigma2;

%sigma2=10^(-10);    %用户引入热噪声
%sigmar2=10^(-10);   %active RIS引入热噪声

f_c=5;      %载频

K=4;            %用户数量
N=512;          %RIS单元数
eta_k=ones(K,1);%权值

large_fading_AI=2.2;
large_fading_DI=2.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsum=zeros(length(user_X),ExNumber);
Rsum_noRIS=zeros(length(user_X),ExNumber);
Rsum_random=zeros(length(user_X),ExNumber);
Rsum_passive=zeros(length(user_X),ExNumber);

for a=1:length(user_X)
    fprintf('第%d轮\n',a);
    parfor b=1:ExNumber
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser]=Position_generate(K,user_X(a));    %基站RIS用户位置设置
        [ h_k,f_k,G] = Channel_generate(K,N,M,large_fading_AI,large_fading_DI,Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser,f_c);

        Theta=diag(exp(1j*2*pi*rand(N,1)));
        W=exp(1j*2*pi*rand(K*M,1))*sqrt(Ps_max/K/M);
        
        [W, Rsum_noRIS(a,b)]= NoRIS_precoding(M,K,N,Ps_max,sigma2,eta_k,W,h_k,f_k,G);
        
%        Rsum_noRIS(a,b)=0;
        [W, ~, Rsum_random(a,b)]= random_RIS_precoding(M,K,N,Ps_max,sigma2,eta_k,Theta,W,h_k,f_k,G);
        [W, Theta, Rsum_passive(a,b)]= passive_RIS_precoding(M,K,N,Ps_max,sigma2,eta_k,Theta,W,h_k,f_k,G);
        Theta=100*Theta;
%        [W,Theta,Rsum(a,b)]= active_RIS_precoding(M,K,N,Ps_max/2,Pr_max/2,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G);
        [W,Theta,Rsum(a,b)]= active_RIS_precoding(M,K,N,Ps_max*0.99,Pr_max*0.01,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G);
        
    end
end

Rsum_mean=mean(Rsum,2);
Rsum_noRIS_mean=mean(Rsum_noRIS,2);
Rsum_passive_mean=mean(Rsum_passive,2);
Rsum_random_mean=mean(Rsum_random,2);

figure;
hold on;
box on;
grid on;
plot(user_X,Rsum_mean,'-r^','LineWidth',1.5);
plot(user_X,Rsum_passive_mean,'-bo','LineWidth',1.5);
plot(user_X,Rsum_random_mean,'-m^','LineWidth',1.5);
plot(user_X,Rsum_noRIS_mean,'--k','LineWidth',1.5);
ylabel('Sum-rate (bps/Hz)','Interpreter','latex');
xlabel('Distance $L$ (m)','Interpreter','latex');

set(gca,'FontName','Times','FontSize',12);

ylim([0 40]);
%legend('Active RIS~~\,~($P_{\rm BS}^{\rm max}=5$ W, $P_{\rm A}^{\rm max}=5$ W)','Passive RIS~\,~($P_{\rm BS}^{\rm max}=10$ W)','Random phase shift ($P_{\rm BS}^{\rm max}=10$ W)','Without RIS ($P_{\rm BS}^{\rm max}=10$ W)','Interpreter','latex','FontSize',12);

legend('Active RIS (ideal case)','Passive RIS [23]','Random phase shift [40]','Without RIS [40]','Interpreter','latex','FontSize',12);

%legend('Active RIS~~\,~($P_{\rm BS}=5$ W, $P_{\rm A}=5$ W)','Passive RIS~\,~($P_{\rm BS}=10$ W)','Without RIS ($P_{\rm BS}=10$ W)','Interpreter','latex','FontSize',12);

save('main.mat','Rsum_mean','Rsum_passive_mean','Rsum_random_mean','Rsum_noRIS_mean');


toc