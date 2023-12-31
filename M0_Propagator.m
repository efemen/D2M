clear; clc; close all;

%% Assumptions

% All planetary orbits are coplanar.
% Earth is a sphere.
% No plane change except from launch to equatorial orbit.


run("M1_Launch.m")

run("M2_Injection.m")

run("M3_Interplanetary.m")

run("M4_MarsCapture.m")
