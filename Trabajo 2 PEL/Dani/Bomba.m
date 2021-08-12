classdef Bomba < handle
    % Objeto bomba
    
    properties (Access = public)
        eta           % Rendimiento
        ro              % Denisdad del liquido
        Pd              % Presion deposito
        pi_perdidas     % Perdidas de presion
        delta_P         % Salto de presion en la bomba
        Ps              % Presion salida bomba
        tau         % Trabajo especifico de la bomba
        Pc
    end
    
    methods
        
        function obj = Bomba(eta, ro, Pd, pi_perdidas) % CONSTRUCTOR
            obj.eta = eta;
            obj.ro = ro;
            obj.Pd = Pd;
            obj.pi_perdidas = pi_perdidas;
            
        end
        
        function [] = Punto_Funcionamiento(obj, Pc)
            obj.Pc = Pc;
            obj.delta_P = (Pc/obj.pi_perdidas) - obj.Pd;
            obj.Ps = obj.delta_P + obj.Pd;
            obj.tau = obj.delta_P./(obj.eta*obj.ro);
        end

    end
end

