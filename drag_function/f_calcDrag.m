function [drag_N,speed_Ma,pDyn_Pa,P_static_Pa,coef_drag]=f_calcDrag(V_ECEF_mag_km_s,altitude_km,ref_area_m2,alpha_deg)
% Calculates drag force on a spacecraft.
% Bahadir Onur Guduru
% 19.11.2023
arguments
V_ECEF_mag_km_s = 0 % Velocity Magnitude ECEF in km/s
altitude_km  = 0 % Altitude from surface in km
ref_area_m2 = 1 % Referance area for the rocket in m^2
alpha_deg   = 0 % Angle of attack in degrees )not 
end

[T,a,P_static_Pa,rho,nu]=atmosisa(altitude_km*1000, "extended", "on", "action", "None"); 
pDyn_Pa=0.5*rho*(V_ECEF_mag_km_s*1000)^2;

% Use drag data from Sutton's book. See "create_drag_data.m"
persistent filedata  
if isempty(filedata)
filedata = load("drag_data.mat");
end

speed_Ma=V_ECEF_mag_km_s*1000/a;

if speed_Ma < 5
    coef_drag=interp2(filedata.drag_data.alpha_deg,filedata.drag_data.Speed_Ma,filedata.drag_data.drag_coef_table,alpha_deg,speed_Ma);
else
    speed_Ma = 5;
    coef_drag=interp2(filedata.drag_data.alpha_deg,filedata.drag_data.Speed_Ma,filedata.drag_data.drag_coef_table,alpha_deg,speed_Ma);
    speed_Ma = V_ECEF_mag_km_s*1000/a;
end

drag_N=ref_area_m2*coef_drag*pDyn_Pa;
end