clear;
load launchdata.mat

%% FP - AOA
figure(1)
plot(T(1:264), fp_angle_ecef(1:264), 'LineWidth', 3)
hold on
xlabel("Time (s)")
ylabel("Flight Path Angle (deg)")
ylim([0 92])

yyaxis right
plot(T(1:264), aoa(1:264), 'LineWidth', 3)
ylabel("Angle of Attack (deg)")
ylim([0 10])
grid on

ax = gca;
ax.FontSize = 12;
set(gcf, 'Position',  [600, 400, 800, 600])

legend("Flight Path Angle", "Angle of Attack")

%% CD - Time
figure(2)
plot(T(1:264), cd_hist(1:264), 'LineWidth', 3)
hold on
xlabel("Time (s)")
ylabel("C_D (-)")
% ylim([0 92])
% 
% yyaxis right
% plot(T(1:264), aoa(1:264), 'LineWidth', 3)
% ylabel("Angle of Attack (deg)")
% ylim([0 10])
grid on

ax = gca;
ax.FontSize = 12;
set(gcf, 'Position',  [600, 400, 800, 600])

% legend("Flight Path Angle", "Angle of Attack")