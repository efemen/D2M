classdef CelestialObject
    % This class defines a celestial object in heliocentric trajectory 
    % and contains the physical properties.
    
    properties
        name
        mass
        r
        r_soi
        r_orbit
        w_orbit
        T_orbit
        mu
        mu_sun
        t_anomaly
    end
    
    methods
        
        function obj = CelestialObject(name, mass, r, r_orbit, jdt)
            % Constants
            m_sun = 1.989e30;  % kg
            G = 6.6743015e-20; %  km^3 k^-1 s^-2
            mu_sun = m_sun * G;

            % Initialize properties
            obj.mass = mass;
            obj.mu = mass * G;
            obj.r = r;
            obj.r_orbit = r_orbit;
            obj.r_soi = obj.r_orbit * (mass / m_sun)^(2/5);
            obj.T_orbit = 2 * pi * obj.r_orbit^1.5 / sqrt(mu_sun);
            obj.w_orbit = sqrt(mu_sun/r_orbit) / r_orbit;
            obj.name = name;
            obj.mu_sun = mu_sun;
            
            
            % True anomaly calculation
            pos = planetEphemeris(jdt, "SolarSystem", obj.name);
            obj.t_anomaly = rad2deg(angle(pos(1)+ 1i * pos(2)));
            
            if obj.t_anomaly < 0
                obj.t_anomaly = obj.t_anomaly + 360;
            end

        end



    end
end

