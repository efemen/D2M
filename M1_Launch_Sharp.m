% Launch will take place from Taiwan. December 2043.

lat = 25.5;
lon = 121.5;

jdt = juliandate([2043, 12, 8, 16, 44, 38]);
sdt = mod(siderealTime(jdt) + lon, 360);

uf = UtilityFunctions();
earth = CelestialObject("Earth", 5.9722e24, 6371.0, 149.598e6, 23.44, jdt);

%% Setup Geometry and Plots
% Earth Sphere
[X,Y,Z] = sphere;

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
n_sun =  uf.ECI2ICRF((-earth.heliocentric_pos / earth.r_orbit)');

% Plots

figure(1);
earth_map = surf(X_E,Y_E,-Z_E);
earthMap = imread("world_Map.jpg");
set(earth_map,'CData', earthMap,'FaceColor','texturemap',"EdgeColor","none")
hold on
colormap white
axis equal
set(gca,'Color','#BEBEBE');
set(gca, 'GridColor', 'white'); 
rotate(earth_map, [0 0 1], siderealTime(jdt))
view(240, 30)
ecliptic_plane = surf(X_ecliptic, Y_ecliptic, Z_ecliptic);
ecliptic_plane.FaceAlpha = 0.2;
hold on

quiver3(6378,0,0, 4000, 0, 0,"filled","LineWidth", 3,"ShowArrowHead","on", "Color","green","MaxHeadSize",10);
text(6378 * 2,0,0,"Vernal Eq. ♈")

quiver3(0, 0, 0, 1e4*n_ecliptic(1), 1e4*n_ecliptic(2), 1e4*n_ecliptic(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "blue");
text(1e4*n_ecliptic(1), 1e4*n_ecliptic(2), 1e4*n_ecliptic(3), "Ecliptic North Pole")

quiver3(0, 0, 0, 1e4*n_sun(1), 1e4*n_sun(2), 1e4*n_sun(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "yellow");
text(1e4*n_sun(1), 1e4*n_sun(2), 1e4*n_sun(3), "Sun ☉")

clear X_E Y_E Z_E X Y Z X_ecliptic Y_ecliptic Z_ecliptic

%% Actual Launch

A = 90;
rho = [sind(A), cosd(A), 0];
V_hat = uf.t_xX(lat, sdt) * rho';

R_i = [earth.r * sind(90 - lat) * cosd(sdt), earth.r * sind(90 - lat) * sind(sdt), earth.r * cosd(90 - lat)];

earth_w = 2 * pi / 8.61640905e4;
magV = 2 * pi * earth.r * cosd(lat) / 8.61640905e4; % km/s
V_i = (magV * V_hat)';


%% Setup Rocket Initial Conditions

dt = 1;          % seconds
T = 0:dt:1200;    % Time matrix
N = length(T);   % Iteration length
m = zeros(N, 1); % Mass matrix
alt = m;         % km
fp_angle = m;    % deg
fp_angle_ecef = m; % deg

alt(1) = 0; % km

X_R = zeros(N, 3);
V_R = X_R;
A_R = V_R;

X_R(1,:) = R_i;
V_R(1,:) = V_i;

% Rocket Parameters -- RD-180, ATLAS 401

m0 = 333e3; % kg
mp = 284e3; % kg


m1 = m0 - mp; %kg
m(1) = m0;    % kg
m2 = m1 - 20830; % kg

g0 = 9.81; % m/s^2
I_sp = 320; % s
c = I_sp * g0; % m/s
m_dot = 1250; % kg/s
T_mag = m_dot * c; % N

% Second Stage Parameters -- Centaur Upper Stage

I_sp2 = 450; % s
c2 = I_sp * g0; % m/s
m_dot2 = 22.453; % kg/s
T_mag2 = m_dot * c; % N
second_stage_flag = 0; % Boolean


roll_program_threshold = 25.5; % km
second_stage_threshold = 500;  % km
p_orbit_speed = sqrt(earth.mu / (earth.r + second_stage_threshold));

a = @(X, V, m)  T_mag * (X / norm(X)) / m / 1e3 - earth.mu * X / norm(X)^3; % Thrust aligned w/ zenith.


%% Iteration


for i = 1:N
    % Timing
    t = T(i);

    % Altitude update
    alt(i) = norm(X_R(i, :)) - earth.r;
    V_ECEF = V_R(i, :) - V_i;
    
    fp_angle(i) = acosd(dot(X_R(i, :), V_R(i, :)) / (norm(X_R(i,:)) * norm(V_R(i,:))));
    fp_angle_ecef(i) = 90 - acosd(dot(X_R(i, :), V_ECEF) / (norm(X_R(i,:)) * norm(V_ECEF)));

    if alt(i) > roll_program_threshold
        a = @(X, V, m)  T_mag * (V / norm(V)) / m / 1e3 - earth.mu * X / norm(X)^3;
        disp("Roll program initated! @ timestep " + string(i))
        roll_program_threshold = 1e24;
    end
    
    % RK4 Solver
    [X_R, V_R] = uf.RK4_launch(a, m(i), dt, X_R, V_R, i);
    A_R(i, :) = (V_R(i + 1) - V_R(1)) / dt;

    % Plot current position.
    figure(1)
    plot3(X_R(i,1), X_R(i, 2), X_R(i, 3),".","Color","#FF3131");
    pause(0.001)
    rotate(earth_map, [0 0 1], rad2deg(earth_w*dt), c_Rot)

    % Staging and burn control.

    if m(i) <= m1
        disp("MECO! @ timestep " + string(i))
        dt = 5;
        disp(string(alt(i)) + " km")
        T_mag = 0;
        a = @(X, V, m)  T_mag * uf.hat(V - V_i) / m / 1e3 - earth.mu * X / norm(X)^3; % T_mag update.
        m_dot = 0;
        m1 = 0;
    end

    if abs(alt(i) - second_stage_threshold) < 1 && norm(V_R(i, :)) < p_orbit_speed*0.5 
        disp("Second stage burn! @ timestep " + string(i))
        disp(string(alt(i)) + " km")
        dt = 1;
        second_stage_flag = 1;
        m_dot = m_dot2;
    end

    if second_stage_flag == 1
        [ra, dec] = uf.ECI2raDec(X_R(i, :));
        burn_direction = (uf.t_xX(dec, ra) * rho')';
        a = @(X, V, m)  T_mag2 * burn_direction / m / 1e3 - earth.mu * X / norm(X)^3; % T_mag update.
    end

    if m(i) <= m2 || norm(V_R(i, :)) >= (p_orbit_speed - 0.1) && second_stage_flag == 1
        second_stage_flag = 0;
        disp("Second stage burnout @ timestep " + string (i))
        dt = 40;
        disp(string(alt(i)) + " km")
        T_mag = 0;
        a = @(X, V, m)  T_mag / m / 1e3 - earth.mu * X / norm(X)^3; % T_mag update.
        m_dot = 0;
        m2 = 0;
    end

    if norm(X_R(i, :)) < earth.r
        disp("Earth impact!")
        break
    end

    % Mass and time update
    m(i + 1) = m(i) - m_dot * dt;
    T(i + 1) = T(i) + dt;
    earth.refresh(t);
    
end

