clear;
load altitude2.mat

%% Gravity
figure(1)

plot(T(1:386), g_hist(1:386)*1e3 .* sind(fp_angle_ecef(1:386)), "LineWidth", 3)



area(T(1:386), g_hist(1:386)* 1e3 .* sind(fp_angle_ecef(1:386)), 'FaceColor', 'black', 'FaceAlpha', 0.1)
ylim([0 10.5])
xlim([0 895])
grid on
xlabel("Time (s)")
ylabel("g sin(\theta) (m/s^2)")
title("Gravitational Acceleration")
ax = gca;
ax.FontSize = 12;
% set(gcf, 'Position',  [600, 400, 1000, 600])
grid on

trapz(T(1:386), g_hist(1:386)*1e3 .* sind(fp_angle_ecef(1:386)))

% figure(2)
% 
% area(T(1:251), drag_hist(1:251) ./ m(1:251), 'FaceColor', 'black', 'FaceAlpha', 0.1)
% xlabel("Time (s)")
% ylabel("Aerodynamic acceleration (m/s^2)")
% grid on
% trapz(T(1:251), drag_hist(1:251) ./ m(1:251))
% 
% yyaxis right
% plot(T(1:251), q_hist(1:251) / 1e3, 'LineWidth', 3)
% ylabel("P_{dyn} (kPa)")