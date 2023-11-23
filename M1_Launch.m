clear; clc; close all;

addpath("drag_function/");

% Launch will take place from East Pacific. 8 December 2043.

lat = 23.43928;
lon = 132.9424;
launch_date = datetime([2043, 12, 8, 15, 59, 0]);
injection_date = datetime([2043, 12, 30, 15, 45, 0]);

jdt = juliandate(launch_date);

sdt = mod(siderealTime(jdt) + lon, 360);


earth = CelestialObject("Earth", 5.9722e24, 6371.0, 149.598e6, 23.43928, jdt);
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
n_ecliptic = uf.t_xX(90-earth.tilt, 270) * [0; 0; 1]; % RA = 270 deg, Dec = (90 - Earth tilt) deg

[X_ecliptic, Y_ecliptic] = meshgrid(-8000:2000:8000); % Generate x and y data
Z_ecliptic = -1/n_ecliptic(3) * (n_ecliptic(1)*X_ecliptic + n_ecliptic(2)*Y_ecliptic); % Solve for z data

% Sun Direction
n_sun =  uf.ICRF2ECI((uf.hat(-earth.heliocentric_pos)'));

% Earth Velocity Direction
n_earth_velocity = uf.rodrigues_rot(n_sun, n_ecliptic, -90);

% Plots

figure(1);
set(gcf, 'Position',  [1200, 0, 800, 800])
earth_map = surf(X_E,Y_E,-Z_E);
earthMap = imread("world_Map.jpg");
set(earth_map,'CData', earthMap,'FaceColor','texturemap',"EdgeColor","none")
hold on
colormap white
axis equal
set(gca, 'Color', 'none')
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
fp_angle_ecef = m; % deg
aoa = m;           % deg
energy = m;        % Specific orbital energy, km^2 / s^2
drag_hist = m;
q_hist = m;
cd_hist = m;
Ma_hist = m;

X_R = zeros(N, 3);
V_R = X_R;
A_R = V_R;

X_R(1, :) = R_i;
V_R(1, :) = V_i;

% % Rocket Parameters -- RD-180, ATLAS 401
% 
% m0 = 333e3; % kg
% mp = 284e3; % kg
% 
% m1 = m0 - mp; %kg
% m1_inert = 21054;
% m(1) = m0;    % kg
% m2 = m1 - m1_inert - 20830; % kg
% 
% g0 = 9.81; % m/s^2
% I_sp = 320; % s
% c = I_sp * g0; % m/s
% m_dot = 1250; % kg/s
% T_mag = m_dot * c; % N
% 
% % Second Stage Parameters -- Centaur Upper Stage
% 
% I_sp2 = 450; % s
% c2 = I_sp * g0; % m/s
% m_dot2 = 22.453; % kg/s
% T_mag2 = m_dot * c2; % N

% LV Parameters -- RD-180, ATLAS 401

m0 = 331781.70; % Launch mass, kg
mp1 = 284089;   % First stage propellant, kg
mp2 = 20830;    % Second stage propellant, kg
m1_inert = 21054; % kg

m1 = m0 - mp1; % MECO mass,
m(1) = m0;    % kg
m2 = m1 - m1_inert - mp2; % SECO mass, kg

g0 = 9.81; % m/s^2
I_sp = 311; % s
c = I_sp * g0; % m/s
m_dot = 1253; % kg/s
T_mag = m_dot * c; % N

% Second Stage Parameters -- Centaur Upper Stage

I_sp2 = 450; % s
c2 = I_sp * g0; % m/s
m_dot2 = 22; % kg/s
T_mag2 = m_dot * c2; % N


second_stage_flag = 0; % Boolean
roll_program_flag = 0; % Boolean

roll_program_threshold = 0.5; % km
burn_direction = uf.hat(X_R(1, :)); % Unit  vector

r_parking_orbit = earth.r + 500; % km
v_parking_orbit = sqrt(earth.mu / r_parking_orbit);
w_parking_orbit = rad2deg(v_parking_orbit / r_parking_orbit);
e_parking_orbit = -0.5 * earth.mu / r_parking_orbit;

a = @(X, V, m, drag)  (T_mag * uf.hat(X) + drag) / m / 1e3 - earth.mu * X / norm(X)^3; % Thrust aligned w/ zenith.

g_sum = 0;
aero_sum = 0;
dv_sum = 0;

%% Iteration

for i = 1:N
    % Timing
    t = T(i);

    % Altitude update
    alt(i) = norm(X_R(i, :)) - earth.r;
    
    [ra_r, dec_r] = uf.ECI2raDec(X_R(i, :));

    V_ECEF = V_R(i, :) - V_i;
    
    fp_angle(i) = 90 - uf.angle_between(X_R(i, :), V_R(i, :));
    fp_angle_ecef(i) = 90 - uf.angle_between(X_R(i, :), V_ECEF);
    aoa(i) = uf.angle_between(burn_direction, V_ECEF);
   

    % Drag Calculation

    [drag, Ma_hist(i), q_hist(i), Ps, cd_hist(i)]  = f_calcDrag(norm(V_ECEF), alt(i), (4 * pi), aoa(i));
    drag  = -drag * uf.hat(V_ECEF);
    drag_hist(i) = norm(drag);

    % Losses
    g_sum = g_sum + norm(earth.mu * X_R(i, :) / norm(X_R(i, :))^3 * dt);
    aero_sum = aero_sum + norm(drag) * dt / m(i) / 1e3;
    dv_sum = dv_sum + T_mag * dt / m(i) / 1e3;
    

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

    if m(i) <= m1 && second_stage_flag == 0
        disp("MECO! @ timestep " + string(i))
        disp(string(alt(i)) + " km")
        a = @(X, V, m, drag) - earth.mu * X / norm(X)^3; % T_mag update.
        m(i) = m(i) - m1_inert;
        m1 = 0;
        aero_loss = aero_sum;
        disp("Second stage burn! @ timestep " + string(i))
        disp(string(alt(i)) + " km")
        dt = 0.5;
        second_stage_flag = 1;
        m_dot = m_dot2;
        T_mag = T_mag2;
    end

    if second_stage_flag == 1
        burn_direction = uf.hat(V_R(i, :));
        a = @(X, V, m, drag)  (T_mag * burn_direction + drag) / m / 1e3 - earth.mu * X / norm(X)^3; % T_mag update.
    end


    if (orbit_now.r_apoapsis - earth.r) > 490 && second_stage_flag == 1
        disp("Coasting until circularization burn... @ timestep " + string (i))
        second_stage_flag = 2;
        dt = 7;
        disp(string(alt(i)) + " km")
        m_dot = 0;
        T_mag = 0;
        a = @(X, V, m, drag) - earth.mu * X / norm(X)^3; % T_mag update.
        g_loss = g_sum;
    end

    if second_stage_flag == 2 && abs(alt(i) - 550) < 3
        disp("Circularization burn! @ timestep = " + num2str(i))
        disp(alt(i));
        second_stage_flag = 3;
        dt = 0.5;
        m_dot = m_dot2;
        T_mag = T_mag2;
        figure(1)
        view(0, 90 - earth.tilt)
    end

    if second_stage_flag == 3
        burn_direction = uf.hat(cross(uf.hat(orbit_now.h_vector), uf.hat(X_R(i, :))));
        a = @(X, V, m, drag)  (T_mag * burn_direction) / m / 1e3 - earth.mu * X / norm(X)^3; % T_mag update.
        

        if abs(energy(i) - e_parking_orbit) < 0.15
            a = @(X, V, m, drag) - earth.mu * X / norm(X)^3; % T_mag update.
            disp("Circularization complete! @ timestep = " + num2str(i))
            dt = 10;
            disp("Perigee height " + string(orbit_now.r_periapsis - earth.r) + " km")
            disp("Apogee height " + string(orbit_now.r_apoapsis - earth.r) + " km")
            second_stage_flag = 4;
            m_dot = 0;
            V_R(i + 1, :) = V_R(i,:);
            T_mag = 0;
        end
    end
    

    if second_stage_flag == 4
        figure(1)
        view(ra_r+sdt, 30);
    end

    if second_stage_flag == 4 && abs(alt(i) - 500) < 1 && abs(norm(V_R(i, :)) - v_parking_orbit) > 0.01
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
        dv_sum = dv_sum + dv_ideal;
        disp("Gravity loss = " + string(g_loss))
        disp("Aerodynamic loss = " + string(aero_loss))
        disp("dV = " + string(dv_sum))
        disp("dV - Gravity loss - Aero loss + V_site = " + string(dv_sum - g_loss - aero_loss + norm(V_i)))
        disp("Parking orbit speed = " + string(v_parking_orbit))
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
% 
% save launchdata.mat

time_to_injection = injection_date - launch_date + seconds(t);  % Time left to injection date in seconds.
dw = mod(w_parking_orbit * seconds(time_to_injection), 360);    % True anomaly change during parking orbit.

X_injection = uf.rodrigues_rot(X_R(i, :), uf.hat(orbit_now.h_vector), dw);
V_injection = uf.rodrigues_rot(V_R(i, :), uf.hat(orbit_now.h_vector), dw);

save parking_orbit.mat X_injection V_injection injection_date dv_sum

