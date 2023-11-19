% Tabulates graph data acquired from image: mach_cd_graph_sutton.jpg by grabit.
% Bahadir Onur Guduru
% 19.11.2023


clear all

% Load
load("drag_graph_data.mat")

% Plot and inspect the graph data
plot(aoa0(:,1),aoa0(:,2))
hold on
plot(aoa4(:,1),aoa4(:,2))
plot(aoa6(:,1),aoa6(:,2))
plot(aoa8(:,1),aoa8(:,2))
plot(aoa10(:,1),aoa10(:,2))
ylabel('Drag Coefficient [CD]')
title('Drag for Rockets')
xlabel('Speed [Ma]')
grid on; grid minor;
hold off

% Tabulate

drag_data.Speed_Ma=linspace(0,6);

drag_data.drag_coef_table=round([interp1(aoa0(:,1),aoa0(:,2),drag_data.Speed_Ma,"linear","extrap")
    interp1(aoa4(:,1),aoa4(:,2),drag_data.Speed_Ma,"linear","extrap")
    interp1(aoa6(:,1),aoa6(:,2),drag_data.Speed_Ma,"linear","extrap")
    interp1(aoa8(:,1),aoa8(:,2),drag_data.Speed_Ma,"linear","extrap")
    interp1(aoa10(:,1),aoa10(:,2),drag_data.Speed_Ma,"linear","extrap")],2)';

drag_data.alpha_deg=[0 4 6 8 10];

% Try
alpha_deg=3
Speed_Ma=1.2
coef_drag=interp2(drag_data.alpha_deg,drag_data.Speed_Ma,drag_data.drag_coef_table,alpha_deg,Speed_Ma)

% Save
save("drag_data.mat","drag_data")


