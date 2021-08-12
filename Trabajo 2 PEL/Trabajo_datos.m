clear all; close all; clc

%% DATOS

% Densidades [kg/m^3]
rho_ox = 1140; 
rho_f = 71;

% Presiones [Pa]
P_s = 40000; % P adaptado
P_amb = 101325; % {SL}
%P_amb = 0; % {vacío}
%P_amb = 40000;

% Presiones de los depósitos [Pa]
P_dox = 6*10^5; 
P_df = 3*10^5; 

% Rendimientos
rend_bf = 0.62;    
rend_bo = 0.64;
rend_t = 0.72;
rend_mec = 0.94;

% Caídas de presión
Pi_iny_cc_ox = 0.78;    
Pi_iny_cc_red = 0.82;
Pi_refr = 0.88;
Pi_iny_pc_ox = 0.74;
Pi_pc = 0.92;
Pi_t = 12;

% Temperaturas [K]
T_ref = 298.15;
T_pc = 900;
T_tin = 900;

% Masas molares [g/mol] 
M_h20 = (2*1.008 + 15.999);
M_o2 = 2*15.999;
M_h2 = 2*1.008;

% Cp [J/molK]
Cp_h20 = 54.32;
Cp_o2 = 39.07;
Cp_h2 = 36.11;

    %Unidades [J/kgK]
    Cp_h20_kg = 54.32/(M_h20*1e-3);
    Cp_o2_kg = 39.07/(M_o2*1e-3); 
    Cp_h2_kg = 36.11/(M_h2*1e-3);



% Constantes
%Entalp_h2o = - 241818; % [J/mol]
%Entalp_h2o_kg = - 241818 * 1000/M_h20; % [J/kg]
Entalp_h2o = -241818e3/(M_h20); %[J/kg]
R_un = 8.3144598; %[J/molK]
g0 = 9.80665; % m/s^2

% Datos Empuje
E_kg = 50*1e3; % kg
E = E_kg*g0; % N

% Datos Pesos
K_o=0.02;
K_f=0.2;
r=0.2;