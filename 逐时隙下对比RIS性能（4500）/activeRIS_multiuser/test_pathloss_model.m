% 测试双射线路损模型
% 验证断点距离和分段特性

clear; clc;

% 参数设置
lambda = 0.06;  % 波长 (m), 对应 5GHz
h_BS = 40;      % 基站高度 (m)
h_RIS = 3;      % RIS高度 (m)
h_User = 8;     % 用户高度 (m)

% 计算断点距离
d_break_BU = (4 * h_BS * h_User) / lambda;
d_break_BR = (4 * h_BS * h_RIS) / lambda;
d_break_RU = (4 * h_RIS * h_User) / lambda;

fprintf('=== 双射线模型断点距离 ===\n');
fprintf('BS-User 断点: %.2f km\n', d_break_BU/1000);
fprintf('BS-RIS 断点:  %.2f km\n', d_break_BR/1000);
fprintf('RIS-User 断点: %.2f km\n', d_break_RU/1000);
fprintf('\n');

% 测试距离范围
d_test = logspace(2, 4.5, 100);  % 100m 到 ~31km

% 计算三条链路的路损
PL_BU = zeros(size(d_test));
PL_BR = zeros(size(d_test));
PL_RU = zeros(size(d_test));

for i = 1:length(d_test)
    PL_BU(i) = getPathLoss_TwoRay(d_test(i), h_BS, h_User, lambda);
    PL_BR(i) = getPathLoss_TwoRay(d_test(i), h_BS, h_RIS, lambda);
    PL_RU(i) = getPathLoss_TwoRay(d_test(i), h_RIS, h_User, lambda);
end

% 转换为dB
PL_BU_dB = 10*log10(PL_BU);
PL_BR_dB = 10*log10(PL_BR);
PL_RU_dB = 10*log10(PL_RU);

% 绘图
figure('Position', [100, 100, 1200, 400]);

% 子图1: BS-User
subplot(1,3,1);
semilogx(d_test/1000, PL_BU_dB, 'b-', 'LineWidth', 2); hold on;
xline(d_break_BU/1000, 'r--', 'LineWidth', 1.5, 'Label', sprintf('断点 %.1f km', d_break_BU/1000));
grid on; box on;
xlabel('距离 (km)', 'FontSize', 11);
ylabel('路损增益 (dB)', 'FontSize', 11);
title('BS-User 直射链路', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.1, 30]);
legend('路损曲线', 'Location', 'southwest', 'FontSize', 10);

% 子图2: BS-RIS
subplot(1,3,2);
semilogx(d_test/1000, PL_BR_dB, 'g-', 'LineWidth', 2); hold on;
xline(d_break_BR/1000, 'r--', 'LineWidth', 1.5, 'Label', sprintf('断点 %.1f km', d_break_BR/1000));
% 标注仿真范围
patch([6, 10, 10, 6]/1, [-200, -200, 0, 0], 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
text(8, -50, '仿真范围', 'FontSize', 10, 'HorizontalAlignment', 'center');
grid on; box on;
xlabel('距离 (km)', 'FontSize', 11);
ylabel('路损增益 (dB)', 'FontSize', 11);
title('BS-RIS 反射前段', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.1, 30]);
legend('路损曲线', 'Location', 'southwest', 'FontSize', 10);

% 子图3: RIS-User
subplot(1,3,3);
semilogx(d_test/1000, PL_RU_dB, 'm-', 'LineWidth', 2); hold on;
xline(d_break_RU/1000, 'r--', 'LineWidth', 1.5, 'Label', sprintf('断点 %.1f km', d_break_RU/1000));
% 标注优化核心区
patch([0.5, 5, 5, 0.5], [-200, -200, 0, 0], 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
text(2, -50, '优化核心区', 'FontSize', 10, 'HorizontalAlignment', 'center');
grid on; box on;
xlabel('距离 (km)', 'FontSize', 11);
ylabel('路损增益 (dB)', 'FontSize', 11);
title('RIS-User 反射后段', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.1, 30]);
legend('路损曲线', 'Location', 'southwest', 'FontSize', 10);

sgtitle('双射线路损模型：分段特性验证', 'FontSize', 14, 'FontWeight', 'bold');

% 验证连续性
fprintf('=== 断点处连续性验证 ===\n');
for link = 1:3
    if link == 1
        h1 = h_BS; h2 = h_User; name = 'BS-User';
    elseif link == 2
        h1 = h_BS; h2 = h_RIS; name = 'BS-RIS';
    else
        h1 = h_RIS; h2 = h_User; name = 'RIS-User';
    end
    
    d_b = (4 * h1 * h2) / lambda;
    PL_before = getPathLoss_TwoRay(d_b * 0.9999, h1, h2, lambda);
    PL_at = getPathLoss_TwoRay(d_b, h1, h2, lambda);
    PL_after = getPathLoss_TwoRay(d_b * 1.0001, h1, h2, lambda);
    
    fprintf('%s: 断点前=%.6e, 断点处=%.6e, 断点后=%.6e\n', ...
        name, PL_before, PL_at, PL_after);
    fprintf('  相对误差: %.4f%%\n', abs(PL_before - PL_after) / PL_at * 100);
end
fprintf('\n');

% 验证衰减指数
fprintf('=== 衰减指数验证 ===\n');
d1 = 500; d2 = 1000;  % 近场测试距离
PL1 = getPathLoss_TwoRay(d1, h_RIS, h_User, lambda);
PL2 = getPathLoss_TwoRay(d2, h_RIS, h_User, lambda);
alpha_near = -log10(PL2/PL1) / log10(d2/d1);
fprintf('RIS-User 近场 (%.0f-%.0f m): 衰减指数 = %.2f (理论值: 2.00)\n', d1, d2, alpha_near);

d1 = 5000; d2 = 10000;  % 远场测试距离
PL1 = getPathLoss_TwoRay(d1, h_BS, h_RIS, lambda);
PL2 = getPathLoss_TwoRay(d2, h_BS, h_RIS, lambda);
alpha_far = -log10(PL2/PL1) / log10(d2/d1);
fprintf('BS-RIS 远场 (%.0f-%.0f m): 衰减指数 = %.2f (理论值: 4.00)\n', d1, d2, alpha_far);

fprintf('\n测试完成！\n');
