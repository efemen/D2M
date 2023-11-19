classdef LaunchVehicle
    % This class defines the launch vehicle wth its mass and all staging
    % details to clean up the code.
    
    properties
        m0 % launch mass
        m_meco % meco mass
        m1 % second stage burn mass
        m_seco % second stage burnout mass
        I_sp1
        I_sp2
        m_dot1
        m_dot2
    end
    
    methods
        function obj = LaunchVehicle(m0, m_meco, m1, m_seco, I_sp2, m_dot1, m_dot2)
            % Initialize the rocket object.
            obj.m0 = m0;
            obj.m1 = m1;
            obj.m_meco = m_meco;
            obj.m_seco = m_seco;
            obj.I_sp1 = I_sp1;
            obj.I_sp2 = I_sp2;
            obj.m_dot1 = m_dot1;
            obj.m_dot2 = m_dot2;

        end
        
    end
end

