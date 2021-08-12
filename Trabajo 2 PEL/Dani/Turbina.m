classdef Turbina < handle
    % Objeto Turbina
    
    properties
        eta       % Rendimiento
        pi        % Relacion expansion
        Pe        % Presion de entrada
        Ps        % Presion de salida
        Te        % Temperatura de entrada
        Ts        % Temperatura de salida

        gamma     % gamma gas       
        cp        % cp gas
        
        tau       % Trabajo especifico de la turbina
        
    end
    
    methods
        function obj = Turbina(eta, pi, Te, gamma, cp) % CONSTRUCTOR
            obj.eta = eta;
            obj.pi = pi;
            % obj.Pe = Pe;
            % obj.Ps = Pe/pi;   
            obj.Te = Te;
            obj.gamma = gamma;
            obj.cp = cp;
            
            obj.Ts = Te*( 1 - eta*( 1 - (1/pi).^((gamma-1)./gamma) ) );
            obj.tau = eta.*cp*Te.*( 1 - (1/pi).^( (gamma-1)./gamma ) );
        end
        
        function [] = Punto_Funcionamiento(obj, Pe)
            obj.Pe = Pe;
            obj.Ps = Pe/obj.pi;
        end
    end
end
