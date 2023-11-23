clear; clc; close all;

%% Assumptions

% All planetary orbits are coplanar.
% mars is a sphere.
% No plane change except from launch to equatorial orbit.

load SOI_IN.mat

%% Object Initilization

uf = UtilityFunctions();

mars = CelestialObject("Mars", 0.64169e24, 3389.5, 227.956e6, 0, jdt_mars); % name, mass, planetary radius, heliocentric radius, jdt


%% Setup Geometry and Plots
% Mars Sphere
[X,Y,Z] = sphere;

X_E = X * mars.r;
Y_E = Y * mars.r;
Z_E = Z * mars.r;
c_RotX = mean(mean(X_E));
c_RotY = mean(mean(Y_E));
c_Rot = [c_RotX c_RotY 0];

% Ecliptic Plane
n_ecliptic = [0; 0; 1]; % RA = 270 deg, DEC = 66.56 deg

[X_ecliptic, Y_ecliptic] = meshgrid(-8000:2000:8000); % Generate x and y data
Z_ecliptic = -1/n_ecliptic(3) * (n_ecliptic(1)*X_ecliptic + n_ecliptic(2)*Y_ecliptic); % Solve for z data

% Sun Direction
n_sun =  uf.ICRF2ECI((uf.hat(-mars.heliocentric_pos)'));

% mars Velocity Direction
n_mars_velocity = uf.rodrigues_rot(n_sun, n_ecliptic, -90);

% Plots

figure(1);
set(gcf, 'Position',  [500, 800, 800, 800])
mars_map = surf(X_E,Y_E,-Z_E);
marsMap = imread("mars_Map.jpg");
set(mars_map,'CData', marsMap,'FaceColor','texturemap',"EdgeColor","none")
hold on
colormap white
axis equal
set(gca,'Color','none');
set(gca, 'GridColor', 'none'); 
set(gca,'Visible','off');


view(240, 30)
% 
% ecliptic_plane = surf(X_ecliptic, Y_ecliptic, Z_ecliptic);
% ecliptic_plane.FaceAlpha = 0.2;
% hold on

quiver3(mars.r, 0,0, 4000, 0, 0,"filled","LineWidth", 3,"ShowArrowHead","on", "Color","green","MaxHeadSize",10);
text(mars.r * 2,0,0,"Vernal Eq. ♈")
% 
% quiver3(0, 0, 0, 1e4*n_ecliptic(1), 1e4*n_ecliptic(2), 1e4*n_ecliptic(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "blue");
% text(1e4*n_ecliptic(1), 1e4*n_ecliptic(2), 1e4*n_ecliptic(3), "Ecliptic North Pole")
% 
% quiver3(0, 0, 0, 1e4*n_sun(1), 1e4*n_sun(2), 1e4*n_sun(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "yellow");
% text(1e4*n_sun(1), 1e4*n_sun(2), 1e4*n_sun(3), "Sun ☉")
% 
% quiver3(0, 0, 0, 1e4*n_mars_velocity(1), 1e4*n_mars_velocity(2), 1e4*n_mars_velocity(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "red");
% text(1e4*n_mars_velocity(1), 1e4*n_mars_velocity(2), 1e4*n_mars_velocity(3), "Mars Velocity")


clear X_E Y_E Z_E X Y Z X_ecliptic Y_ecliptic Z_ecliptic



%% Setup Spacecraft Initial Conditions

dt = 500;          % seconds
T = 0:dt:2e6;     % Time matrix
N = length(T);    % Iteration length


X_SC = zeros(N, 3);
V_SC = X_SC;
A_SC = V_SC;


X_i = X_SC_mars;
V_i = V_SC_mars;

mars_w = 2 * pi / 88642.663;

capture_flag = 0;

X_SC(1,:) = X_i;
V_SC(1,:) = V_i;

e = zeros(N,1);
u = e;
ke = u;

a = @(X) -mars.mu * X / norm(X)^3;

for i = 1:776
    A_SC(i,:) = a(X_SC(i,:));
    
    T(i + 1) = T(i) + dt;
    
    [X_SC, V_SC] = uf.RK4(a, dt, X_SC, V_SC, i);
    
    
    % KE, PE and SE calculated as a metric of accuracy.
    u(i) = -mars.mu / norm(X_SC(i,:));
    ke(i) = 0.5 * norm(V_SC(i,:))^2;
    e(i) = u(i) + ke(i);
    orbit_now = OrbitalElements(X_SC(i, :), V_SC(i, :), mars.mu);
    
    
    if norm(X_SC(i,:)) > mars.r_soi + 1e6
      disp("SOI radius reached.")
      soi_timestep = i;
      break
    end

    if abs(norm(X_SC(i,:)) - orbit_now.r_periapsis) < 100 && capture_flag == 0        
        xlim([-2e4 2e4])
        ylim([-2e4 2e4])
        disp("Periapsis reached.")
        V_ideal = uf.hat(uf.rodrigues_rot(X_SC(i + 1, :), [0, 0, 1], 90)) * sqrt(mars.mu / norm(X_SC(i, :)));
        dV = V_ideal - V_SC(i, :);
        V_SC(i, :) =  V_ideal;
        V_SC(i + 1, :)  = V_SC(i, :);
        disp("Capture burn complete. dV = " + string(norm(dV)) + " km/s")
        disp("Circular orbit at r = " + string(norm(X_SC(i, :)))+ " km")
        dv_sum = dv_sum + norm(dV);
        disp("Total mission dV = " + string(dv_sum) + " km/s")

        capture_flag = 1;
    end

   % Plot current position.
    figure(1)
    plot3(X_SC(i,1), X_SC(i, 2), X_SC(i, 3), ".","Color","#FF3131");
    rotate(mars_map, [0 0 1], rad2deg(mars_w*dt), c_Rot)


    if norm(X_SC(i, :)) < mars.r - 1
        disp("Mars impact!")
        break
    end

end

% T_SOI = T(soi_timestep);
% V_SOI = V_SC(soi_timestep, :);
% X_SOI = X_SC(soi_timestep, :);
% disp("Time elapsed since launch is " + string(T_SOI) + " seconds.")
% 
% X_SOI = uf.ECI2ICRF(X_SOI');
% V_SOI = uf.ECI2ICRF(V_SOI');
% jdt_SOI =  juliandate(datetime(injection_date) + seconds(T_SOI));


