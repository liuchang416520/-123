% clear; clc; close all;
% % ==========================================================
% % 45km 鍦烘櫙锛氫笌 main_1000 绠楁硶閫昏緫涓€鑷达紙BCD+SCA 杞ㄨ抗浼樺寲锛夛紝缁撴灉瓒嬪娍璧板悜鐩歌繎
% % 鍙傛暟鎸?45km 灏哄害鍚堢悊璁剧疆锛氳矾寰勬崯鑰楀ぇ闇€ N=256锛屾椂闅欐暟鍦ㄧ簿搴︿笌鑰楁椂闂存姌涓?% % ==========================================================
% 
% % 1. 杞ㄨ抗璧风粓鐐癸紙浣嶄簬鍩虹珯涓嬫柟锛屽嚑浣曠粨鏋勪笌 main_1000 绫讳技锛氬簳閮ㄧ洿绾库啋鍚戠敤鎴峰集鏇诧級
% x_start = 5000;  y_start = -8000;
% x_end   = 40000; y_end   = -8000;
% 
% % 2. 杈圭晫闄愬埗
% y_min = y_start; y_max = 0;      % y 鈭?[y_start, 0]锛屽悜鐢ㄦ埛鏂瑰悜鎺㈢储
% 
% % 3. 杩愬姩瀛﹀弬鏁?% % 鍗曟椂闅欎綅绉?v_max*dt 闇€瑕嗙洊璺緞锛氭按骞崇害 35km + 鍨傜洿绾?8km锛屾€昏矾寰?~43km
% % Nslot=100锛氭瘡鏃堕殭 ~430m锛屽湪杞ㄨ抗绮惧害涓庤繍琛屾椂闂撮棿鎶樹腑锛堝師 150 鏃堕殭鍗曟杩唬绾?5h锛?% v_max = 25;      % 鏈€澶ч€熷害 (m/s)
% Nslot   = 100;   % 鏃堕殭鏁帮細杩囧皯杞ㄨ抗绮楃硻锛岃繃澶氳€楁椂鍓у
% T_total = 2000;  % 鎬绘椂闂?(s)锛宒t=20s锛寁_max*dt=500m
% dt      = T_total / Nslot;
% 
% fprintf('=== 45km 鍦烘櫙 (绠楁硶涓?main_1000 涓€鑷达紝鍙傛暟鎸夊昂搴﹀悎鐞嗚缃? ===\n');
% fprintf('璧风偣: [%.0f, %.0f], 缁堢偣: [%.0f, %.0f]\n', x_start, y_start, x_end, y_end);
% fprintf('Nslot=%d, 鍗曟椂闅欎綅绉? %.0f m\n', Nslot, v_max * dt);
% 
% % 5. 鐢ㄦ埛浣嶇疆锛堜腑蹇冨湪 [45000, 0, 10]锛屽崐寰?1500m锛?% user_X_center = 45000;
% user_Y_center = 0;
% user_R        = 1500;
% 
% % 6. 楂樺害璁剧疆
% ht_BS   = 35;          % 鍩虹珯楂樺害 35m
% h_RIS   = 15;          % USV (RIS) 楂樺害 15m
% hr_User = 10;          % 鐢ㄦ埛楂樺害 10m
% 
% % 7. 閫氫俊鍙傛暟
% % 45km 璺緞鎹熻€楄緝 1km 楂樼害 20*log10(45)鈮?3dB锛岄渶鏇村 RIS 鍗曞厓琛ュ伩
% % main_1000 鐢?N=64 鍥犲叾涓?1km 寰満鏅紱45km 鍦烘櫙 N=256 涓哄悎鐞嗛厤缃?% K = 4; M = 10; N = 256;
% P_total_dBm = 47;       % 澶ц窛绂婚渶楂樺姛鐜?% P_total_mW  = 10^(P_total_dBm/10);
% RIS_power_ratio = 0.1; % 涓?main_1000 涓€鑷?% Pr_max = P_total_mW * RIS_power_ratio;      
% Ps_max = P_total_mW * (1 - RIS_power_ratio); 
% fprintf('鍔熺巼鍒嗛厤: 鍩虹珯 %.2f W, RIS %.2f W\n', Ps_max/1000, Pr_max/1000);
% 
% % 鍣０涓庝俊閬撳弬鏁?% sigma2  = 1e-13; sigmar2 = 1e-13; eta_k = ones(K,1);   
% f_c = 5; lambda = 3e8/(f_c*1e9); K_rician = 10;         
% tau_power = 4; kappa_power = 4; gamma_reflect = 0;
% scale_G = 1e5; scale_f = 1e8;
% 
% % 8. 鍒濆杞ㄨ抗鐢熸垚锛堢嚎鎬ф彃鍊硷級
% x_traj = linspace(x_start, x_end, Nslot).';
% y_traj = linspace(y_start, y_end, Nslot).';
% 
% bs_pos = [0, 0, ht_BS];
% 
% % 鐢ㄦ埛鐢熸垚
% rng(42);
% theta = 2*pi*rand(K,1); rad = user_R*rand(K,1);
% user_pos = zeros(K,3);
% for k = 1:K
%     user_pos(k,:) = [user_X_center + rad(k)*cos(theta(k)), user_Y_center + rad(k)*sin(theta(k)), hr_User];
% end
% 
% % 淇￠亾鍒濆鍖?% rng(123);
% h_k0 = zeros(M,K); f_k0 = zeros(N,K); G0 = zeros(N,M);
% Theta0 = diag(exp(1j*2*pi*rand(N,1)));
% W0 = exp(1j*2*pi*rand(K*M,1))*sqrt(Ps_max/K/M);
% 
% % 缁撴瀯浣撴墦鍖呴厤缃?% traj_cfg = struct();
% traj_cfg.enable = true;
% traj_cfg.bs_pos = bs_pos;
% traj_cfg.user_pos = user_pos;
% traj_cfg.x_traj = x_traj;
% traj_cfg.y_traj = y_traj;
% traj_cfg.h_const = h_RIS;
% traj_cfg.Nslot = Nslot;
% traj_cfg.dt = dt;
% traj_cfg.v_max = v_max;
% traj_cfg.y_min = y_min; traj_cfg.y_max = y_max; % 浼犲叆 y_max=0
% traj_cfg.lambda = lambda;
% traj_cfg.gamma_reflect = gamma_reflect;
% traj_cfg.K_rician = K_rician;
% traj_cfg.tau_power = tau_power;
% traj_cfg.kappa_power = kappa_power;
% traj_cfg.no_direct_link = false;
% traj_cfg.MaxIter = 15;  % 涓?main_1000 涓€鑷?% traj_cfg.MinIter = 5;
% traj_cfg.allow_hover_at_max_wsr = false;
% traj_cfg.return_to_endpoint = false;
% traj_cfg.user_X_center = user_X_center; traj_cfg.user_R = user_R; 
% traj_cfg.user_Y_center = user_Y_center; 
% traj_cfg.scale_G = scale_G; traj_cfg.scale_f = scale_f;
% traj_cfg.use_smart_init = false;   % 澶у満鏅笅鍏抽棴锛岄伩鍏?warm start 澶辨晥鏃?E1/E2/E3 涓嶇ǔ瀹?% traj_cfg.trust_step = 0.5 * v_max * dt;  % SCA 淇¤禆鍩熶笌鍗曟椂闅欎綅绉诲尮閰嶏紝閬垮厤杩囧ぇ姝ラ暱
% 
% % ==========================================================
% % 璋冪敤鏍稿績浼樺寲鍑芥暟
% % ==========================================================
% [W_opt, Theta_opt, Rsum, ~, ~, x_opt, y_opt, rate_hist, ~] = active_RIS_precoding(M, K, N, Ps_max, Pr_max, sigma2, sigmar2, eta_k, Theta0, W0, h_k0, f_k0, G0, traj_cfg);
% 
% % ==========================================================
% % 缁撴灉缁樺浘
% % ==========================================================
% wsr_final = rate_hist(:, end);
% [wsr_max, best_slot] = max(wsr_final);
% best_x = x_opt(best_slot); best_y = y_opt(best_slot);
% high_wsr_idx = find(wsr_final >= 0.8 * wsr_max);
% 
% figure('Position', [100, 100, 1200, 500]);
% 
% % 瀛愬浘1锛氳建杩瑰浘锛堜笌 main_1000 椋庢牸涓€鑷达級
% subplot(1, 2, 1); hold on; grid on; box on;
% plot(bs_pos(1), bs_pos(2), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
% scatter(user_pos(:,1), user_pos(:,2), 50, 'b', 'filled');
% th = linspace(0, 2*pi, 100);
% plot(user_X_center + user_R*cos(th), user_Y_center + user_R*sin(th), '--', 'Color', [0.6 0.6 0.6]);
% 
% plot(x_traj, y_traj, '--k', 'LineWidth', 1.5);  % 鍒濆杞ㄨ抗 (鐩寸嚎)
% plot(x_opt, y_opt, '-r', 'LineWidth', 2, 'Marker', '.');
% % 鏃堕殭鏍囨敞锛堟瘡闅?5 涓椂闅欙紝涓?main_1000 涓€鑷达級
% step_mark = 5;
% for n = 1:step_mark:Nslot
%     y_offset = 500;
%     text(x_opt(n), y_opt(n) + y_offset, sprintf('t=%d', n), 'FontSize', 9, 'Color', [0.5 0 0.5], 'FontWeight', 'bold');
%     plot(x_opt(n), y_opt(n), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [0.5 0 0.5]);
% end
% if mod(Nslot-1, step_mark) ~= 0
%     text(x_opt(end), y_opt(end) + 500, sprintf('t=%d', Nslot), 'FontSize', 9, 'Color', [0.5 0 0.5], 'FontWeight', 'bold');
% end
% 
% plot(x_start, y_start, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
% plot(x_end, y_end, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% 
% legend({'BS', 'Users', 'User Region', '鍒濆杞ㄨ抗 (鐩寸嚎)', '浼樺寲鍚庤建杩?}, 'Location', 'best');
% xlabel('x (m)'); ylabel('y (m)');
% title('USV杞ㄨ抗鍥?(45km 鍦烘櫙绠楁硶椹卞姩鍔涙祴璇?', 'FontSize', 13, 'FontWeight', 'bold');
% xlim([-1000, 46000]); 
% ylim([-9000, 2000]); 
% yline(0, 'k-', 'y=0 杈圭晫', 'LineWidth', 1.5); 
% 
% % 瀛愬浘2锛歐SR锛堜笌 main_1000 涓€鑷达紝鏍囧嚭宄板€硷級
% subplot(1, 2, 2); hold on; grid on; box on;
% plot(1:Nslot, wsr_final, '-r', 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
% plot(best_slot, wsr_max, 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
% text(best_slot, wsr_max + 0.3, sprintf('Max: %.2f', wsr_max), 'HorizontalAlignment', 'center');
% xlabel('鏃堕殭鏁?); ylabel('WSR (bps/Hz)');
% title(sprintf('鍔犳潈鍜岄€熺巼鍙樺寲 (骞冲潎 WSR: %.4f)', mean(wsr_final)), 'FontSize', 13, 'FontWeight', 'bold');
% 
% fprintf('\n45km 鍦烘櫙锛氳建杩瑰簲鍚?y=0 (鐢ㄦ埛渚? 寮洸锛學SR 鍦ㄩ潬杩戠敤鎴峰杈惧埌宄板€笺€俓n');


clear; clc; close all;
% ==========================================================
% 鍦烘櫙瀹氫箟锛氫笅鏂圭洿绾垮垏鍏?+ 閫愭椂闅?Slot-by-slot)璐績瀵讳紭
% 鎵╁睍鑷?45km 澶у昂搴﹀満鏅?% ==========================================================

% 1. 杞ㄨ抗璧风偣
x_start = 35000;  y_start = -5000;

% 2. 杈圭晫闄愬埗锛氫粎淇濈暀 y 鏂瑰悜绾︽潫
y_min = -6000;   y_max = 0;      % y 鈭?[y_min, 0]

% 3. 杩愬姩瀛﹀弬鏁?v_max = 25;      % 鏈€澶ч€熷害 (m/s)
T_total = 5500;  % 鎬婚琛屾椂闂?5000绉?(绾?3鍒嗛挓)
Nslot   = 200;   % 鏃堕殭鏁?150
dt      = T_total / Nslot;

fprintf('=== 閫愭椂闅欏浼樺弬鏁?===\n');
fprintf('璧风偣: [%.0f, %.0f]\n', x_start, y_start);
fprintf('鎬绘椂闂? %d s, 鏃堕殭鏁? %d, 鍗曟椂闅? %.2f s\n', T_total, Nslot, dt);
fprintf('鏈€澶у崟姝ユ満鍔ㄨ窛绂? %.1f m\n', v_max * dt);

% 4. 鐢ㄦ埛浣嶇疆锛堜腑蹇冨湪 [45000, 0, 10]锛屽崐寰?1500m锛?user_X_center = 45000;
user_Y_center = 0;
user_R        = 1500;

% 5. 楂樺害璁剧疆
ht_BS   = 35;          % 鍩虹珯楂樺害 35m
h_RIS   = 15;          % USV (RIS) 楂樺害 15m
hr_User = 10;          % 鐢ㄦ埛楂樺害 10m

% 6. 閫氫俊鍙傛暟
K = 4; M = 10; N = 256;
P_total_dBm = 47;
P_total_mW  = 10^(P_total_dBm/10);
RIS_power_ratio = 0.05; % 鍔熺巼姣斾緥
Pr_max = P_total_mW * RIS_power_ratio;      
Ps_max = P_total_mW * (1 - RIS_power_ratio); 
fprintf('鍔熺巼鍒嗛厤: 鍩虹珯 %.2f W, RIS %.2f W\n', Ps_max/1000, Pr_max/1000);

% 鍣０涓庝俊閬撳弬鏁?sigma2  = 1e-13; sigmar2 = 1e-13; eta_k = ones(K,1);   
f_c = 5; lambda = 3e8/(f_c*1e9); K_rician = 10;         
tau_power = 4; kappa_power = 4; gamma_reflect = 0;

% 銆愬叧閿慨鏀广€戯細璋冨皬缂╂斁鍊嶆暟锛岄槻姝?45km 绾ц仈淇￠亾瀵艰嚧 E绯绘暟 鐖嗙偢(10^22)
scale_G = 1e2; scale_f = 1e4;

% 7. 鍒濆杞ㄨ抗鐢熸垚锛堛€愬叧閿慨鏀广€戯細鍏ㄩ儴鍒濆鍖栦负璧风偣锛?x_traj = x_start * ones(Nslot, 1);
y_traj = y_start * ones(Nslot, 1);

bs_pos = [0, 0, ht_BS];

% 鐢ㄦ埛鐢熸垚
rng(42);
theta = 2*pi*rand(K,1); rad = user_R*rand(K,1);
user_pos = zeros(K,3);
for k = 1:K
    user_pos(k,:) = [user_X_center + rad(k)*cos(theta(k)), user_Y_center + rad(k)*sin(theta(k)), hr_User];
end

% 淇￠亾鍒濆鍖?rng(123);
h_k0 = zeros(M,K); f_k0 = zeros(N,K); G0 = zeros(N,M);
Theta0 = diag(exp(1j*2*pi*rand(N,1)));
W0 = exp(1j*2*pi*rand(K*M,1))*sqrt(Ps_max/K/M);

% 缁撴瀯浣撴墦鍖呴厤缃?traj_cfg = struct();
traj_cfg.enable = true;
traj_cfg.bs_pos = bs_pos;
traj_cfg.user_pos = user_pos;
traj_cfg.x_traj = x_traj;
traj_cfg.y_traj = y_traj;
traj_cfg.h_const = h_RIS;
traj_cfg.Nslot = Nslot;
traj_cfg.dt = dt;
traj_cfg.v_max = v_max;
traj_cfg.y_min = y_min; traj_cfg.y_max = y_max; 
traj_cfg.lambda = lambda;
traj_cfg.gamma_reflect = gamma_reflect;
traj_cfg.K_rician = K_rician;
traj_cfg.tau_power = tau_power;
traj_cfg.kappa_power = kappa_power;
traj_cfg.no_direct_link = false;
traj_cfg.MaxIter = 18;  
traj_cfg.MinIter = 10;
traj_cfg.disp_thresh = 500;
% 銆愬叧閿慨鏀广€戯細鍏抽棴寮哄埗杩斿洖缁堢偣绾︽潫鍜屽己琛屾偓鍋滃垽瀹?traj_cfg.allow_hover_at_max_wsr = false;
traj_cfg.return_to_endpoint = false;
traj_cfg.user_X_center = user_X_center; traj_cfg.user_R = user_R; 
traj_cfg.user_Y_center = user_Y_center; 
traj_cfg.scale_G = scale_G; traj_cfg.scale_f = scale_f;
traj_cfg.use_smart_init = false;
traj_cfg.e2_ratio_cap = 2.5;
traj_cfg.e_clip = 4.0;

% ==========================================================
% 璋冪敤閫愭椂闅欐牳蹇冧紭鍖栧嚱鏁?% ==========================================================
[W_opt, Theta_opt, Rsum, ~, ~, x_opt, y_opt, rate_hist, ~] = active_RIS_precoding(M, K, N, Ps_max, Pr_max, sigma2, sigmar2, eta_k, Theta0, W0, h_k0, f_k0, G0, traj_cfg);

% ==========================================================
% 缁撴灉缁樺浘
% ==========================================================
wsr_final = rate_hist; % 鐜板湪鏄崟缁村害鏁扮粍
[wsr_max, best_slot] = max(wsr_final);

figure('Position', [100, 100, 1200, 500]);

% 瀛愬浘1锛氳建杩瑰浘
subplot(1, 2, 1); hold on; grid on; box on;

% 鐢诲熀绔欎笌鐢ㄦ埛
plot(bs_pos(1), bs_pos(2), 'ks', 'MarkerSize', 14, 'MarkerFaceColor', [0.3 0.8 0.3], 'LineWidth', 2);
scatter(user_pos(:,1), user_pos(:,2), 80, 'b', 'filled', 'MarkerEdgeColor', 'k');
th = linspace(0, 2*pi, 200);
plot(user_X_center + user_R*cos(th), user_Y_center + user_R*sin(th), '--', 'Color', [0.6 0.6 0.6]);

% 鐢讳紭鍖栧悗鐨勮建杩?plot(x_opt, y_opt, '-r', 'LineWidth', 2.5, 'Marker', '.', 'MarkerSize', 10);

% 鍔ㄦ€佹墦鏃堕殭鏍囩 (姣?10 涓椂闅欐爣涓€娆?
step_mark = 10;
for n = 1:step_mark:Nslot
    text(x_opt(n), y_opt(n) + 300, sprintf('t=%d', n), 'FontSize', 9, 'Color', [0.5 0 0.5]);
    plot(x_opt(n), y_opt(n), 'o', 'MarkerSize', 5, 'MarkerEdgeColor', [0.5 0 0.5]);
end

% 鏍囨敞璧风偣
plot(x_start, y_start, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
text(x_start-1500, y_start, '璧风偣', 'FontSize', 11, 'FontWeight', 'bold');

xlabel('x (m)'); ylabel('y (m)');
title('USV 閫愭椂闅欒椽蹇冨浼樿建杩?, 'FontSize', 13, 'FontWeight', 'bold');
xlim([-1000, 48000]); 
ylim([-7000, 2000]); 
yline(0, 'k-', 'y=0 杈圭晫', 'LineWidth', 1.5); 
legend({'BS', 'Users', 'User Region', 'Trajectory', 'Start'}, 'Location', 'best');

% 瀛愬浘2锛歐SR
subplot(1, 2, 2); hold on; grid on; box on;
plot(1:Nslot, wsr_final, '-r', 'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
plot(best_slot, wsr_max, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y');
text(best_slot, wsr_max + max(wsr_final)*0.05, sprintf('Max: %.2f (t=%d)', wsr_max, best_slot), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

xlabel('鏃堕殭'); ylabel('WSR (bps/Hz)');
title(sprintf('閫愭椂闅欏姞鏉冨拰閫熺巼 (Peak: %.4f)', wsr_max), 'FontSize', 13, 'FontWeight', 'bold');

fprintf('\n=== 娴嬭瘯瀹屾垚 ===\n');
fprintf('瑙傚療绾㈢嚎鏄惁鍏嬫湇浜嗗熀绔欏紩鍔涳紝鎴愬姛椋炲悜鍙充晶 45km 鐨勭敤鎴峰尯鍩熴€俓n');
