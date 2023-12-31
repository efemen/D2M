classdef LaunchVehicle
    % This class defines the launch vehicle wth its mass and all staging
    % details to clean up the code.
    
    properties
        m     % current mass
        m_dot % current m_dot
        T     % current thrust
        stage % current stage

        m0    % launch mass
        m_p1
        m_p2
        m_inert1 % first stage inert mass
        m_inert2 % second stage inert mass
   
        
        T_mag1 % first stage thrust
        T_mag2 % second stage thrust
        m_dot1 % first stage mass flow rate
        m_dot2 % second stage mass flow rate
        
    end
    
    methods
        function obj = LaunchVehicle(m0, m_p1, m_p2, m_inert1, m_inert2, T_mag1, T_mag2, m_dot1, m_dot2)
            % Initialize the launch vehicle object.
            obj.m = m0;
            obj.m_dot = m_dot1;
            obj.stage = 1;
            obj.burning = 1;

            obj.m0 = m0;
            obj.m_p1 = m_p1;
            obj.m_p2 = m_p2;
            obj.m_inert1 = m_inert1;
            obj.m_inert2 = m_inert2;

            obj.T_mag1 = T_mag1;
            obj.T_mag2 = T_mag2;

            obj.m_dot1 = m_dot1;
            obj.m_dot2 = m_dot2;
        end

        function obj = stop_burn(obj)
            obj.T = 0;
            obj.m_dot = 0;
        end

        function obj = start_burn(obj)
            if obj.stage == 1
                obj.T = obj.T_mag1;
                obj.m_dot = obj.T_mag2;
            else
                obj.T = obj.T_mag2;
                obj.m_dot = obj.m_dot2;
            end
        end


        function obj = stage_sep(obj)
            if obj.stage == 1
                obj.m = obj.m - obj.m_inert1;
                obj.T = obj.T_mag2;
                obj.m_dot = obj.m_dot2;
            else
                obj.m = obj.m - obj.m_inert2;
            end

            obj.stage = obj.stage + 1;
        end
        
        function obj = refresh(obj, dt)
            obj.m = obj.m - obj.m_dot * dt;
        end

    end
end