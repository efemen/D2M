clear; clc;

load SOI_OUT.mat X_SOI T_SOI V_SOI

jdt = juliandate(datetime([2043, 12, 28, 16, 44, 0]) + seconds(T_SOI));

%% Object Initilization

uf = UtilityFunctions();

earth = CelestialObject("Earth", 5.9722e24, 6371.0, 149.598e6, 23.44, jdt); % name, mass, planetary radius, heliocentric radius, jdt
mars = CelestialObject("Mars", 0.64169e24, 3389.5, 227.956e6, 0, jdt); % name, mass, planetary radius, heliocentric radius, jdt

% Spacecraft SOI-Exit conditions
V_SOI(3) = 0;
X_SOI(3) = 0;

X_i =  X_SOI' + earth.heliocentric_pos;
V_i =  V_SOI' + earth.heliocentric_vel;

%% Setup Geometry and Plots
figure(1)

plot(earth.heliocentric_pos(1), earth.heliocentric_pos(2), "b.", "MarkerSize", 10)
hold on
grid on
axis equal
xlim([-3e8 3e8])
ylim([-3e8 3e8])
quiver(earth.heliocentric_pos(1), earth.heliocentric_pos(2), earth.heliocentric_vel(1), earth.heliocentric_vel(2), 1e6, "blue")

plot(mars.heliocentric_pos(1), mars.heliocentric_pos(2), "r.", "MarkerSize", 10)
quiver(mars.heliocentric_pos(1), mars.heliocentric_pos(2), mars.heliocentric_vel(1), mars.heliocentric_vel(2), 1e6, "red")
sun = plot(0, 0, "y.", "MarkerSize", 30);


plot(X_i(1), X_i(2), "r.", "MarkerSize", 10)
quiver(X_i(1), X_i(2), V_i(1), V_i(2), 1e6, "red")


%% Analytical Orbital Elements

V_ri = dot(uf.hat(X_i), V_i);

h_vector = cross(X_i, V_i);
h = norm(h_vector);
e_vector = (1/earth.mu_sun) * ((norm(V_i)^2 - earth.mu_sun / norm(X_i)) * X_i - norm(X_i) * V_ri * V_i);
eccentricity = norm(e_vector);

arg_periapsis = uf.angle_between(e_vector, [1, 0, 0]);

t_anomaly = 0:0.5:360;
r = (h^2 / earth.mu_sun) ./ (1 + eccentricity * cosd(t_anomaly));

X_SC = r' .* [cosd(t_anomaly'), ....
        sind(t_anomaly')];

for i = 1:length(X_SC)
    X_SC(i, :) = [cosd(arg_periapsis), -sind(arg_periapsis); sind(arg_periapsis) cosd(arg_periapsis)]...
                 * X_SC(i, :)';
end


plot(X_SC(:, 1), X_SC(:, 2))
