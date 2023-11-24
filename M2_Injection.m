clear; clc; close all;
addpath("images/")

%% Assumptions

% All planetary orbits are coplanar.
% Earth is a sphere.
% No plane change except from launch to equatorial orbit.


%% User Parameters

load parking_orbit.mat

figure(99)
three_weeks = imshow(imread("3weeks.png"));
set(gcf, 'Position',  [0, 0, 1920, 1080])
pause(3)
close all

jdt = juliandate(injection_date);
sdt = siderealTime(jdt); % Earth map update

%% Object Initilization

uf = UtilityFunctions();

earth = CelestialObject("Earth", 5.97217e24, 6371.0084, 1.49598e8, 23.43928, jdt); % name, mass, planetary radius, heliocentric radius, jdt
mars = CelestialObject("Mars", 0.641691e24, 3389.50, 2.27939e8, 0, jdt); % name, mass, planetary radius, heliocentric radius, jdt


%% Setup Geometry and Plots
% Earth Sphere
[X,Y,Z] = sphere(100);

X_E = X * earth.r;
Y_E = Y * earth.r;
Z_E = Z * earth.r;
c_RotX = mean(mean(X_E));
c_RotY = mean(mean(Y_E));
c_Rot = [c_RotX c_RotY 0];

% Ecliptic Plane
n_ecliptic = uf.t_xX(90-earth.tilt, 270) * [0; 0; 1]; % RA = 270 deg, DEC = 66.56 deg

[X_ecliptic, Y_ecliptic] = meshgrid(-8000:2000:8000); % Generate x and y data
Z_ecliptic = -1/n_ecliptic(3) * (n_ecliptic(1)*X_ecliptic + n_ecliptic(2)*Y_ecliptic); % Solve for z data

% Sun Direction
n_sun =  uf.ICRF2ECI((uf.hat(-earth.heliocentric_pos)'));

% Earth Velocity Direction
n_earth_velocity = uf.rodrigues_rot(n_sun, n_ecliptic, -90);

% Plots

figure(1);
earth_map = surf(X_E,Y_E,-Z_E);
earthMap = imread("world_Map.jpg");
set(earth_map,'CData', earthMap,'FaceColor','texturemap',"EdgeColor","none")
hold on
colormap white
axis equal
rotate(earth_map, [0 0 1], sdt)
uf.draw_space();

ecliptic_plane = surf(X_ecliptic, Y_ecliptic, Z_ecliptic);
ecliptic_plane.FaceAlpha = 0.2;
ecliptic_plane.EdgeColor = "white";
ecliptic_plane.EdgeAlpha = 0.2;
hold on

quiver3(6378,0,0, 4000, 0, 0,"filled","LineWidth", 3,"ShowArrowHead","on", "Color","green","MaxHeadSize",10);
text(6378 * 2,0,0,"Vernal Eq. ♈", "Color", "white")

quiver3(0, 0, 0, 1e4*n_ecliptic(1), 1e4*n_ecliptic(2), 1e4*n_ecliptic(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "blue");
text(1e4*n_ecliptic(1), 1e4*n_ecliptic(2), 1e4*n_ecliptic(3), "Ecliptic North Pole", "Color", "white")

quiver3(0, 0, 0, 1e4*n_sun(1), 1e4*n_sun(2), 1e4*n_sun(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "yellow");
text(1e4*n_sun(1), 1e4*n_sun(2), 1e4*n_sun(3), "Sun ☉", "Color", "white")

quiver3(0, 0, 0, 1e4*n_earth_velocity(1), 1e4*n_earth_velocity(2), 1e4*n_earth_velocity(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "red");
text(1e4*n_earth_velocity(1), 1e4*n_earth_velocity(2), 1e4*n_earth_velocity(3), "Earth Velocity 🜨", "Color", "white")
set(gcf, 'Position',  [0, 0, 1920, 1080])

clear X_E Y_E Z_E X Y Z X_ecliptic Y_ecliptic Z_ecliptic



%% Maneuvering Parameters

departureStatus = 0;
SC_R = earth.r + 500;

transfer_e = (mars.r_orbit - earth.r_orbit) / (mars.r_orbit + earth.r_orbit);
perihelion_velocity = sqrt(earth.mu_sun * (1 + transfer_e) / earth.r_orbit);
hyperbolic_excess = perihelion_velocity - earth.w_orbit * earth.r_orbit;
injection_velocity = sqrt(hyperbolic_excess^2 + 2 * earth.mu / SC_R);
injection_dv = injection_velocity - sqrt(earth.mu / SC_R);
departure_angle = 180 + acosd(1 / (1 + (SC_R * hyperbolic_excess^2 / earth.mu)));

departureThreshold = 500; % km
departure_coordinates = SC_R * uf.rodrigues_rot(n_earth_velocity, n_ecliptic, departure_angle)';

timestepThreshold = 2e4; % km
timestepStatus = 0;
extendedTimeStep = 500; % s

%% Setup Spacecraft Initial Conditions

dt = 10;          % seconds
T = 0:dt:2e5;     % Time matrix
N = length(T);    % Iteration length


X_SC = zeros(N, 3);
V_SC = X_SC;
A_SC = V_SC;


earth_w = 7.292115e-5; % rad / s

X_SC(1,:) = X_injection;
V_SC(1,:) = V_injection;
telemetry = text(X_SC(1, 1) + 30, X_SC(1, 2) + 30, X_SC(1,3) + 30, "Alt = " + 500 + " km" + newline + ...
    "V = " + norm(V_SC(1, :))+ " km/s");
telemetry.Color = "white";

e = zeros(N,1);
u = e;
ke = u;
alt = u;

a = @(X) -earth.mu * X / norm(X)^3;

plot3(departure_coordinates(1), departure_coordinates(2), departure_coordinates(3), "o", 'MarkerSize', 20, "MarkerFaceColor", "white");
text(departure_coordinates(1)*1.1, departure_coordinates(2)*1.1, departure_coordinates(3)*1.1, "Injection Burn Point", "Color", "white")

for i = 1:N
    A_SC(i,:) = a(X_SC(i,:));
    
    T(i + 1) = T(i) + dt;
    
    [X_SC, V_SC] = uf.RK4(a, dt, X_SC, V_SC, i);
    
    
    % KE, PE and SE calculated as a metric of accuracy.
    u(i) = -earth.mu / norm(X_SC(i,:));
    ke(i) = 0.5 * norm(V_SC(i,:))^2;
    e(i) = u(i) + ke(i);
    alt(i) = norm(X_SC(i, :)) - earth.r;
    [ra_r, dec_r] = uf.ECI2raDec(X_SC(i, :));

  
    if abs(norm(X_SC(i,:) - departure_coordinates)) < departureThreshold
       disp("Departure time!")
       X_SC(i,:) = departure_coordinates;
       X_SC(i + 1, :) = X_SC(i, :);

       V_SC(i + 1, :) = injection_velocity  * uf.hat(V_SC(i, :));
       V_SC(i, :) =  V_SC(i + 1, :);

       departureThreshold = 0;
       departureStatus = 1;
       disp("Departure trajectory achieved!")
       disp("dV = " + injection_dv + " km/s")
       dv_sum = dv_sum + injection_dv;
    end
    
    if norm(X_SC(i,:)) > timestepThreshold && timestepStatus == 0
       dt = extendedTimeStep;
       disp("Time step has been increased to " + dt + " seconds.")
       ts_left = length(T) - i;
       T(i:end) = T(i):dt:(T(i) + dt*ts_left);
       timestepStatus = 1;
    end

    if norm(X_SC(i, :)) > timestepThreshold * 10 && timestepStatus == 1
        dt = 1000;
    end
    
    if norm(X_SC(i,:)) > earth.r_soi
      disp("SOI radius reached.")
      soi_timestep = i;      
      break
    end


   % Plot current position.
    figure(1)
    plot3(X_SC(i,1), X_SC(i, 2), X_SC(i, 3), ".","Color","#FF3131");
    rotate(earth_map, [0 0 1], rad2deg(earth_w*dt), c_Rot)
    view(ra_r + sdt + 30, 40);
    telemetry.String = "Alt = " + alt(i) + " km" + newline + ...
    "V = " + norm(V_SC(i, :))+ " km/s";
    telemetry.Position = X_SC(i ,:) + [30 ,30 ,30];

end

T_SOI = T(soi_timestep);
V_SOI = V_SC(soi_timestep, :);
X_SOI = X_SC(soi_timestep, :);
disp("Time elapsed since injection date is " + string(T_SOI) + " seconds.")

X_SOI = uf.ECI2ICRF(X_SOI');
V_SOI = uf.ECI2ICRF(V_SOI');
jdt_SOI =  juliandate(datetime(injection_date) + seconds(T_SOI));

save SOI_OUT.mat X_SOI V_SOI T_SOI jdt_SOI dv_sum
