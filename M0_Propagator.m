clear; clc; close all;

%% Assumptions

% All planetary orbits are coplanar.
% Earth is a sphere.
% No plane change except from launch to equatorial orbit.


%% User Parameters

jdt = juliandate([2043, 12, 8, 15, 59, 0]);

launch_lat = 23.5; % Arbitrary for now
launch_lon = 121.5;

sdt = siderealTime(jdt) + launch_lon;

%% Object Initilization

uf = UtilityFunctions();

earth = CelestialObject("Earth", 5.9722e24, 6371.0, 149.598e6, jdt); % name, mass, planetary radius, heliocentric radius, jdt
mars = CelestialObject("Mars", 0.64169e24, 3389.5, 227.956e6, jdt); % name, mass, planetary radius, heliocentric radius, jdt

run('M1_LaunchParameters.m')   % Get launch parameters



%% Earth Ephemeris
[earth_Pos, earth_V] = planetEphemeris(jdt, "SolarSystem", "Earth");

earth_V_hat = earth_V / norm(earth_V);
earth_R_hat = -earth_Pos / norm(earth_Pos);



%% Mars Ephemeris
[mars_Pos, mars_V] = planetEphemeris(jdt, "SolarSystem", "Mars");
mars_Orbit_R = norm(mars_Pos);  % this will be accepted as the circular radius of heliocentric Mars orbit.
mars_R_hat = mars_Pos / norm(mars_Pos);


%% Maneuvering Parameters
planeThreshold = 50;            % desired plane is accepted as reached if SC is closer than the treshold.
planeChangeStatus = 0;          % plane change status
analytical_dV_planeChange = 2 * magV * sind(lat / 2); % analytical value of the plane change delta-v.

departureStatus = 0;
SC_R = 6378 + 500;

transfer_e = (mars.r_orbit - earth.r_orbit) / (mars_Orbit_R + earth.r_orbit);
perihelion_velocity = sqrt(earth.mu_sun * (1 + transfer_e) / earth.r_orbit);
hyperbolic_excess = perihelion_velocity - earth.w_orbit * earth.r_orbit;
injection_velocity = sqrt(hyperbolic_excess^2 + 2 * earth.mu / SC_R);
injection_dv = injection_velocity - sqrt(earth.mu / SC_R);
departure_angle = 180 + acosd(1 / (1 + (SC_R * hyperbolic_excess^2 / earth.mu)));

departureThreshold = 500; % km
departure_coordinates = SC_R * [cosd(departure_angle) -sind(departure_angle); sind(departure_angle) cosd(departure_angle)] * earth_V_hat(1, 1:2)';
departure_coordinates = [departure_coordinates' 0];
departure_coordinates = departure_coordinates / (norm(departure_coordinates) / SC_R);

timestepThreshold = 2e4; % km
timestepStatus = 0;
extendedTimeStep = 500; % s

%% Setup Spacecraft Initial Conditions

dt = 25;          % seconds
T = 0:dt:2e5;     % Time matrix
N = length(T);    % Iteration length


X_SC = zeros(N, 3);
V_SC = X_SC;
A_SC = V_SC;


X_SC(1,:) = R_i;
V_SC(1,:) = V_i;

e = zeros(N,1);
u = e;
ke = u;

a = @(X) -earth.mu * X / norm(X)^3;

for i = 1:N
   A_SC(i,:) = a(X_SC(i,:));

   T(i + 1) = T(i) + dt;

   [X_SC, V_SC] = uf.RK4(a, dt, X_SC, V_SC, i);


   % KE, PE and SE calculated as a metric of accuracy.
   u(i) = -earth.mu / norm(X_SC(i,:));
   ke(i) = 0.5 * norm(V_SC(i,:))^2;
   e(i) = u(i) + ke(i);


   % this if block checks whether plane change conditions are met.
   if abs(X_SC(i,3)) < planeThreshold
       disp("Plane change time!")
       X_SC(i,3) = 0;
       X_SC(i + 1, :) = X_SC(i, :);
       [ra, dec] = uf.ECI2raDec(X_SC(i,:));
       V_hat = uf.t_xX(dec, ra) * [1; 0; 0];
       
       V_SC(i + 1, :) = (V_hat * magV)';
       dV_planeChange = V_SC(i + 1,:) - V_SC(i,:);
       V_SC(i,:) =  V_SC(i + 1, :);
       planeThreshold = 0;
       planeChangeStatus = 1;
       disp("Plane change successful!")
       disp("Achieved dV = " + string(norm(dV_planeChange)) + " km/s")
       disp("Analytical value dV = " + string(analytical_dV_planeChange) + " km/s")
       disp("Difference = " + string(abs(norm(dV_planeChange) - analytical_dV_planeChange)) + " km/s");
   end

   if planeChangeStatus == 1
       if abs(norm(X_SC(i,:) - departure_coordinates)) < departureThreshold
           disp("Departure time!")
           X_SC(i,:) = departure_coordinates;
           X_SC(i + 1, :) = X_SC(i, :);
           [ra, dec] = uf.ECI2raDec(X_SC(i,:));
           V_hat = uf.t_xX(dec, ra) * [1; 0; 0];
           V_SC(i + 1, :) = (V_hat * injection_velocity)';
           V_SC(i,:) =  V_SC(i + 1, :);
           departureThreshold = 0;
           departureStatus = 1;
           disp("Departure trajectory achieved!")
           disp("dV = " + injection_dv + " km/s")
       end

      if norm(X_SC(i,:)) > timestepThreshold && timestepStatus == 0
           dt = extendedTimeStep;
           disp("Time step has been increased to " + dt + " seconds.")
           ts_left = length(T) - i;
           T(i:end) = T(i):dt:(T(i) + dt*ts_left);
           timestepStatus = 1;
      end

      if norm(X_SC(i,:)) > earth.r_soi
          disp("SOI radius reached.")
          soi_timestep = i;
          break
      end
   end

   % Real-time energy plots for monitoring.
   figure(2)
   subplot(3,1,1)
   plot(T(i),e(i),".r");
   xlabel("Time (s)")
   ylabel("SE (km^2/s^2)")
   ylim([-30 5])
   grid on
   hold on
   
   subplot(3,1,2)
   plot(T(i),norm(X_SC(i,:)),'.r')
%    ylim([6700, 7000])
   xlabel("Time (s)")
   ylabel("R (km)")
   grid on
   hold on

   subplot(3,1,3)
   plot(T(i),norm(A_SC(i,:)),'.r')
   xlabel("Time (s)")
   ylabel("Acceleration (km/s^2)")
%    ylim([0,12])
   grid on
   hold on
end

T_SOI = T(soi_timestep);
V_SOI = V_SC(soi_timestep, :);
X_SOI = X_SC(soi_timestep, :);
disp("Time elapsed since launch is " + string(T_SOI) + " seconds.")

