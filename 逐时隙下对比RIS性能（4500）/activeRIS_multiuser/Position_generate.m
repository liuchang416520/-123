function [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser]=Position_generate(K,user_X,ht_BS,h_RIS,hr_User,RIS_X_Center,user_region_radius)
Dis_BStoRIS=0;
Dis_BStoUser=zeros(K,1);
Dis_RIStoUser=zeros(K,1);

% 海上场景位置：BS在海岸线固定点，RIS更靠近BS
BS_position=[0 0 ht_BS];
RIS_position=[RIS_X_Center 0 h_RIS];

user_position=zeros(K,3);
% 船舶用户位于以 (user_X, 0) 为中心、半径 user_region_radius 的圆内，Z=hr_User
user_position(1,:)=[user_X+user_region_radius*cos(2*pi*rand()) user_region_radius*sin(2*pi*rand()) hr_User];
user_position(2,:)=[user_X+user_region_radius*cos(2*pi*rand()) user_region_radius*sin(2*pi*rand()) hr_User];
user_position(3,:)=[user_X+user_region_radius*cos(2*pi*rand()) user_region_radius*sin(2*pi*rand()) hr_User];
user_position(4,:)=[user_X+user_region_radius*cos(2*pi*rand()) user_region_radius*sin(2*pi*rand()) hr_User];
% user_position(5,:)=[dist+4*(rand()-0.5) 4*(rand()-0.5)];
% user_position(6,:)=[dist+4*(rand()-0.5) 4*(rand()-0.5)];

Dis_BStoRIS=distance(BS_position,RIS_position);

for k=1:K
    user_position_temp=reshape(user_position(k,:),3,1);
    Dis_BStoUser(k)=distance(BS_position,user_position_temp);
end

for k=1:K
    user_position_temp=reshape(user_position(k,:),3,1);
    Dis_RIStoUser(k)=distance(RIS_position,user_position_temp);
end

end