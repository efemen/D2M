classdef UtilityFunctions
    % This file organizes the functions that are utilized in this work in a
    % document.
    
    methods
        function obj = UtilityFunctions()
            disp("UtilityFunctions initialized.")
        end
        
        function xX = t_xX(obj, lat, sid)
            % Topocentric to geocentric equatorial
            xX = zeros(3,3);
            xX(1,:) = [-sind(sid) -sind(lat)*cosd(sid) cosd(lat)*cosd(sid)];
            xX(2,:) = [cosd(sid) -sind(lat)*sind(sid) cosd(lat)*sind(sid)];
            xX(3,:) = [0 cosd(lat) sind(lat)];
        end

        function [ra, dec] = ECI2raDec(obj, R)
            % Convert ECI coordinates to right ascension and declination.
            R_hat = R/norm(R);
            dec = asind(R_hat(3));
            if R(2) > 0
                ra = acosd(R_hat(1) / cosd(dec));
            else
                ra = 360 - acosd((R_hat(1) / cosd(dec)));
            end
        end

        
        function [X_next, V_next] = RK4(obj, dydx, dt, X_SC, V_SC, i)
        % RK4 numerical solver
            dv1 = dt * dydx(X_SC(i,:));
            dx1 = dt * V_SC(i,:);
        
            dv2 = dt * dydx(X_SC(i,:) + 0.5 * dx1);
            dx2 = dt * (V_SC(i,:) + 0.5 * dv1);
        
            dv3 = dt * dydx(X_SC(i,:) + 0.5  * dx2);
            dx3 = dt * (V_SC(i,:) +  0.5 * dv2);
        
            dv4 = dt * dydx(X_SC(i,:) + dx3);
            dx4 = dt * (V_SC(i,:) + dv3);
        
            V_SC(i + 1,:) = V_SC(i,:) +  (dv1 + 2*dv2 + 2*dv3 + dv4) / 6;
            X_SC(i + 1,:) = X_SC(i,:) + (dx1 + 2*dx2 + 2*dx3 + dx4) / 6;
        
            X_next = X_SC;
            V_next = V_SC;
        end
    
        function [X_next, V_next] = RK4_launch(obj, dydx, m, dt, X_SC, V_SC, i)
        % RK4 numerical solver
            dv1 = dt * dydx(X_SC(i,:), V_SC(i, :), m);
            dx1 = dt * V_SC(i,:);
        
            dv2 = dt * dydx(X_SC(i,:) + 0.5 * dx1, V_SC(i, :), m);
            dx2 = dt * (V_SC(i,:) + 0.5 * dv1);
        
            dv3 = dt * dydx(X_SC(i,:) + 0.5  * dx2, V_SC(i, :), m);
            dx3 = dt * (V_SC(i,:) +  0.5 * dv2);
        
            dv4 = dt * dydx(X_SC(i,:) + dx3, V_SC(i, :), m);
            dx4 = dt * (V_SC(i,:) + dv3);
        
            V_SC(i + 1,:) = V_SC(i,:) +  (dv1 + 2*dv2 + 2*dv3 + dv4) / 6;
            X_SC(i + 1,:) = X_SC(i,:) + (dx1 + 2*dx2 + 2*dx3 + dx4) / 6;
        
            X_next = X_SC;
            V_next = V_SC;
        end

    end
end

