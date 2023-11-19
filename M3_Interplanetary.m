clear; clc; close all;

load SOI_OUT.mat X_SOI T_SOI V_SOI jdt_SOI

% jdt = juliandate(datetime([2043, 12, 28, 16, 44, 0]) + seconds(T_SOI));
jdt = jdt_SOI;

%% Object Initilization

uf = UtilityFunctions();

earth = CelestialObject("Earth", 5.9722e24, 6371.0, 149.598e6, 23.44, jdt); % name, mass, planetary radius, heliocentric radius, jdt
mars = CelestialObject("Mars", 0.64169e24, 3389.5, 227.956e6, 0, jdt); % name, mass, planetary radius, heliocentric radius, jdt

% Spacecraft SOI-Exit conditions
V_SOI(3) = 0;
X_SOI(3) = 0;

X_i =  X_SOI' + earth.heliocentric_pos;
V_i =  V_SOI' + earth.heliocentric_vel;
V_i = 0.996 * V_i;

%% Setup Geometry and Plots
figure(1)

earth_plot = plot(earth.heliocentric_pos(1), earth.heliocentric_pos(2), "b.", "MarkerSize", 10);
hold on
grid on
axis equal
xlim([-3e8 3e8])
ylim([-3e8 3e8])

quiver(earth.heliocentric_pos(1), earth.heliocentric_pos(2), earth.heliocentric_vel(1), earth.heliocentric_vel(2), 1e6, "blue");

mars_plot = plot(mars.heliocentric_pos(1), mars.heliocentric_pos(2), "r.", "MarkerSize", 10);
quiver(mars.heliocentric_pos(1), mars.heliocentric_pos(2), mars.heliocentric_vel(1), mars.heliocentric_vel(2), 1e6, "red")
sun = plot(0, 0, "y.", "MarkerSize", 30);


SC_plot = plot(X_i(1), X_i(2), "r.", "MarkerSize", 10);
SC_quiver = quiver(X_i(1), X_i(2), V_i(1), V_i(2), 1e6, "red");
set(gcf, 'Position',  [600, 400, 800, 800])

%% Analytical Orbital Elements

orbit_now = OrbitalElements(X_i, V_i, earth.mu_sun);

t_anomaly_now = uf.angle_between(orbit_now.e_vector, X_i);
E = 2 * atan(sqrt((1-orbit_now.e)/(+1+orbit_now.e)) * tand(t_anomaly_now/2));
Me = E - orbit_now.e * sin(E);
t_periapsis = Me * orbit_now.period / (2*pi);

arg_periapsis = uf.angle_between(orbit_now.e_vector, [1, 0, 0]);
peri_rot = [cosd(arg_periapsis), -sind(arg_periapsis); sind(arg_periapsis) cosd(arg_periapsis)];

t_anomaly = 0:0.5:360;

r = (orbit_now.h^2 / earth.mu_sun) ./ (1 + orbit_now.e * cosd(t_anomaly));

X_transfer = r' .* [cosd(t_anomaly'), ....
        sind(t_anomaly')];

X_Mars = mars.r_orbit .* [cosd(t_anomaly'), ....
        sind(t_anomaly')];

X_Earth = earth.r_orbit .* [cosd(t_anomaly'), ....
        sind(t_anomaly')];


for i = 1:length(X_transfer)
    X_transfer(i, :) = peri_rot * X_transfer(i, :)';
end

plot(X_transfer(:, 1), X_transfer(:, 2))
plot(X_Mars(:, 1), X_Mars(:, 2), "k")
plot(X_Earth(:, 1), X_Earth(:, 2), "k")


N = 600;
T = zeros(1, N);
X_SC = zeros(N, 2);
V_SC = X_SC;

dt = 86400/2; % seconds
min_dist = 1e5;

for i = 1:N
    t = T(i);
    Me = 2 * pi * (t + t_periapsis) / orbit_now.period;
    t_fun = @(E) E - orbit_now.e * sin(E) - Me;
    E = fzero(t_fun, Me);
    t_anomaly_now = rad2deg(2 * atan(tan(E/2) * sqrt((1+orbit_now.e)/(1-orbit_now.e))));
    r_now = (orbit_now.h^2 / earth.mu_sun) ./ (1 + orbit_now.e * cosd(t_anomaly_now));
    X_SC(i, :) = peri_rot * r_now * [cosd(t_anomaly_now); sind(t_anomaly_now)];


    SC_plot.XData = X_SC(i, 1); 
    SC_plot.YData = X_SC(i, 2);

    earth = earth.refresh(dt);
    mars = mars.refresh(dt);

    earth_plot.XData = earth.heliocentric_pos(1);
    earth_plot.YData = earth.heliocentric_pos(2);

    mars_plot.XData = mars.heliocentric_pos(1);
    mars_plot.YData = mars.heliocentric_pos(2);

    pause(0.001)

    sc2mars =  mars.heliocentric_pos(1:2) - X_SC(i, :);
    
    if norm(sc2mars) < mars.r_soi
        disp("SOI REACHED")
        dt = 10;
        min_dist = norm(sc2mars);
    end

    T(i + 1) = t + dt;
end





