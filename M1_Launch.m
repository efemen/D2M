clear; clc; close all;

addpath("drag_function/");

% Launch will take place from Taiwan. December 2043.

lat = 23.44;
lon = 132.9424;

jdt = juliandate([2043, 12, 8, 15, 59, 0]);

sdt = mod(siderealTime(jdt) + lon, 360);


earth = CelestialObject("Earth", 5.9722e24, 6371.0, 149.598e6, 23.44, jdt);
uf = UtilityFunctions();

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
n_sun =  uf.ICRF2ECI((uf.hat(-earth.heliocentric_pos)'));

% Earth Velocity Direction
n_earth_velocity = uf.rodrigues_rot(n_sun, n_ecliptic, 90);

% Plots

figure(1);
set(gcf, 'Position',  [1200, 0, 800, 800])
earth_map = surf(X_E,Y_E,-Z_E);
earthMap = imread("world_Map.jpg");
set(earth_map,'CData', earthMap,'FaceColor','texturemap',"EdgeColor","none")
hold on
colormap white
axis equal
set(gca,'Color','#BEBEBE');
set(gca, 'GridColor', 'white'); 
rotate(earth_map, [0 0 1], siderealTime(jdt))

view(185.53, -37.5487);
campos([-1397, 18386,-6889])
camva(9.07)
camtarget([-194.8, 5976, 2695])

% view(240, 30)
ecliptic_plane = surf(X_ecliptic, Y_ecliptic, Z_ecliptic);
ecliptic_plane.FaceAlpha = 0.2;
hold on

quiver3(6378,0,0, 4000, 0, 0,"filled","LineWidth", 3,"ShowArrowHead","on", "Color","green","MaxHeadSize",10);
text(6378 * 2,0,0,"Vernal Eq. ♈")

quiver3(0, 0, 0, 1e4*n_ecliptic(1), 1e4*n_ecliptic(2), 1e4*n_ecliptic(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "blue");
text(1e4*n_ecliptic(1), 1e4*n_ecliptic(2), 1e4*n_ecliptic(3), "Ecliptic North Pole")

quiver3(0, 0, 0, 1e4*n_sun(1), 1e4*n_sun(2), 1e4*n_sun(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "yellow");
text(1e4*n_sun(1), 1e4*n_sun(2), 1e4*n_sun(3), "Sun ☉")

quiver3(0, 0, 0, 1e4*n_earth_velocity(1), 1e4*n_earth_velocity(2), 1e4*n_earth_velocity(3),"filled","LineWidth", 3,"ShowArrowHead", "on", "Color", "red");
text(1e4*n_earth_velocity(1), 1e4*n_earth_velocity(2), 1e4*n_earth_velocity(3), "Earth Velocity 🜨")

figure(2)
set(gcf, 'Position',  [100, 600, 800, 200])

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

dt = 0.5;          % seconds
T = 0:dt:800;      % time matrix
N = length(T);     % iteration length
m = zeros(N, 1);   % Mass matrix
alt = m;           % km
fp_angle = m;      % deg
fp_angle_wind = m; % deg
aoa = m;           % deg
energy = m;        % Specific orbital energy, km^2 / s^2

X_R = zeros(N, 3);
V_R = X_R;
A_R = V_R;

X_R(1, :) = R_i;
V_R(1, :) = V_i;

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
roll_program_flag = 0; % Boolean

roll_program_threshold = 0.5; % km
burn_direction = uf.hat(X_R(1, :)); % Unit  vector

second_stage_threshold = 490;  % km
r_parking_orbit = earth.r + 500; % km
v_parking_orbit = sqrt(earth.mu / r_parking_orbit);
e_parking_orbit = -0.5 * earth.mu / r_parking_orbit;

a = @(X, V, m, drag)  (T_mag * uf.hat(X) + drag) / m / 1e3 - earth.mu * X / norm(X)^3; % Thrust aligned w/ zenith.

%% Iteration

for i = 1:N
    % Timing
    t = T(i);

    % Altitude update
    alt(i) = norm(X_R(i, :)) - earth.r;
    
    [ra_r, dec_r] = uf.ECI2raDec(X_R(i, :));

    V_ECEF = V_R(i, :) - V_i;
    
    fp_angle(i) = 90 - uf.angle_between(X_R(i, :), V_R(i, :));
    fp_angle_wind(i) = 90 - uf.angle_between(X_R(i, :), V_ECEF);
    aoa(i) = uf.angle_between(burn_direction, V_ECEF);

    % Drag Calculation

    drag = -uf.hat(V_ECEF) * f_calcDrag(norm(V_ECEF), alt(i), (4 * pi), aoa(i));

    % Analytical Orbit
    orbit_now = OrbitalElements(X_R(i, :), V_R(i, :), earth.mu);

    energy(i) = 0.5 * norm(V_R(i,:)) ^ 2 - earth.mu / norm(X_R(i, :));
    
    % RK4 Solver
    [X_R, V_R] = uf.RK4_launch(a, m(i), dt, X_R, V_R, drag, i);
    A_R(i, :) = (V_R(i + 1, :) - V_R(i, :)) / dt;

    % Plot current position.
    figure(1)
    plot3(X_R(i,1), X_R(i, 2), X_R(i, 3), ".","Color","#FF3131");
    rotate(earth_map, [0 0 1], rad2deg(earth_w*dt), c_Rot)

    velocity_vector = quiver3(X_R(i, 1), X_R(i, 2), X_R(i, 3), V_R(i, 1), V_R(i, 2), V_R(i, 3), 6e2, "Color", "blue");
    attitude_vector = quiver3(X_R(i, 1), X_R(i, 2), X_R(i, 3), burn_direction(1), burn_direction(2), burn_direction(3), 6e2, "Color", "green");
    wind_frame_vector = quiver3(X_R(i, 1), X_R(i, 2), X_R(i, 3), V_ECEF(1), V_ECEF(2), V_ECEF(3), 6e2, "Color", "magenta");

    pause(0.001)


    %% Staging and ascent control.

    if alt(i) > roll_program_threshold && alt(i) < 50 && roll_program_flag == 0
        disp("Roll program initated! @ timestep " + string(i))
        dt = 1;
        roll_program_flag = 1;
    end

    if roll_program_flag == 1
        burn_direction = uf.hat(uf.rodrigues_rot(uf.hat(V_ECEF), uf.hat(cross(V_ECEF, V_i)), 0.5 * dt));
        a = @(X, V, m, drag)  (T_mag * burn_direction + drag) / m / 1e3 - earth.mu * X / norm(X)^3;

        if uf.angle_between(V_ECEF, V_R(i, :)) <= 1
            disp("Roll program finished @ timestep " + string(i))
            disp(string(alt(i)) + " km")
            disp(string("Flight path angle = " + fp_angle(i)) + " deg")
            roll_program_flag = 2;
            a = @(X, V, m, drag)  (T_mag * uf.hat(V) + drag) / m / 1e3 - earth.mu * X / norm(X)^3;
        end
    end

    if m(i) <= m1
        disp("MECO! @ timestep " + string(i))
        dt = 2;
        disp(string(alt(i)) + " km")
        T_mag = 0;
        a = @(X, V, m, drag) - earth.mu * X / norm(X)^3; % T_mag update.
        m_dot = 0;
        m1 = 0;
        second_stage_threshold = alt(i);
    end

    if abs(alt(i) - second_stage_threshold) < 1 && second_stage_flag == 0
        disp("Second stage burn! @ timestep " + string(i))
        disp(string(alt(i)) + " km")
        dt = 0.5;
        second_stage_flag = 1;
        m_dot = m_dot2;
    end

    if second_stage_flag == 1
        burn_direction = uf.hat(V_R(i, :));
        a = @(X, V, m, drag)  (T_mag2 * burn_direction + drag) / m / 1e3 - earth.mu * X / norm(X)^3; % T_mag update.
    end


    if (orbit_now.r_apoapsis - earth.r) > 490 && second_stage_flag == 1
        disp("Coasting until circularization burn... @ timestep " + string (i))
        figure(1)
        view(240, 30)
        disp(m(i) - m2)
        second_stage_flag = 2;
        dt = 10;
        disp(string(alt(i)) + " km")
        m_dot = 0;
        a = @(X, V, m, drag) - earth.mu * X / norm(X)^3; % T_mag update.
    end

    if second_stage_flag == 2 && abs(alt(i) - 500) < 1
        disp("Circularization burn!")
        second_stage_flag = 3;
        dt = 1;
        m_dot = m_dot2;
    end

    if second_stage_flag == 3
        burn_direction = uf.hat(cross(uf.hat(orbit_now.h_vector), uf.hat(X_R(i, :))));
        a = @(X, V, m, drag)  (T_mag2 * burn_direction) / m / 1e3 - earth.mu * X / norm(X)^3; % T_mag update.
        

        if abs(energy(i) - e_parking_orbit) < 0.4
            a = @(X, V, m, drag) - earth.mu * X / norm(X)^3; % T_mag update.
            disp("Circularization complete!")
            dt = 20;
            disp("Perigee height " + string(orbit_now.r_periapsis - earth.r) + " km")
            disp("Apogee height " + string(orbit_now.r_apoapsis - earth.r) + " km")
            second_stage_flag = 4;
            m_dot = 0;
        end
    end

    if t > 4000 && t < 5000 && abs(alt(i) - 500) < 1 && abs(norm(V_R(i, :)) - v_parking_orbit) > 0.01
        X_R(i, :) = r_parking_orbit * uf.hat(X_R(i, :)); 
        X_R(i + 1, :) = X_R(i, :);
        V_ideal = v_parking_orbit * uf.hat(cross(orbit_now.h_vector, X_R(i, :)));
        dv_ideal = norm(v_parking_orbit * uf.hat(cross(orbit_now.h_vector, X_R(i, :))) - V_R(i, :));
        V_R(i, :) = V_ideal;
        V_R(i + 1, :) = V_R(i, :);
        
        disp("Orbit idealized dv = " + string(dv_ideal*1e3) + " m/s")
        
        orbit_now = OrbitalElements(X_R(i, :), V_R(i, :), earth.mu);
        X_R = X_R(1:i, :);
        V_R = V_R(1:i, :);
        T = T(1:i);
        m = m(1:i);
        alt = alt(1:i);
        
        figure(1)
        view(120, 30)
        break
    end

    if norm(X_R(i, :)) < earth.r - 1
        disp("Earth impact!")
        break
    end

    % Mass and time update
    m(i + 1) = m(i) - m_dot * dt;
    T(i + 1) = T(i) + dt;
    earth = earth.refresh(dt);
    delete(attitude_vector)
    delete(wind_frame_vector)
    delete(velocity_vector)

    figure(2)
    plot(t, alt(i), "r.")
    xlabel("Time (s)")
    ylabel("Altitude (km)")
    ylim([0 1000])
    xlim([0 1000])
    grid on
    hold on
    
end


