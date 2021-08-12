% Ene ste script metemos todos los datos necesarios para el Probl2, he
% querido separar los datos y la resolución para mayor claridad

% Masas molares [kg/mol]
M_H2O = (2*1.008 + 15.999)/1000; 
M_O2 = 2*15.999/1000;
M_H2 = 2*1.008/1000;

% Cp [J/molK]
Cp_H2O = 54.32;
Cp_O2 = 39.07;
Cp_H2 = 36.11;

% Densidades [kg/m^3]
densidad_ox = 1140;
densidad_f = 71;
densidad_H2O = 1000;

% Constantes [J/molK]
Ru = 8314.4598*1e-3; % J/molK
g0 = 9.80665; % m/s^2


T_ref = 298.15; % K
T_max_t = 900; % K

OF_st = M_O2/2/M_H2; % OF en masas

Entalp_H2O = -241818; % J/mol

% Rendimientos
rend_b_f = 0.62;
rend_b_ox = 0.64;
rend_t = 0.72;
rend_mec = 0.94;

% Caídas de presión
pi_iny_cc_ox = 0.78;
pi_iny_cc_f = 0.82;
pi_ref = 0.88;
pi_iny_pc_ox = 0.74;
pi_pc = 0.92;
pi_t_max = 1/12;

% Presiones de los depósitos
P_dep_ox = 6*1e5; % Pa
P_dep_f = 3*1e5; % Pa

% Datos de empuje:
E_kg = 50*1e3; % kg
E = E_kg*g0; % kg*m/s^2 = N
P_s = 40*1e3; % Pa

% Datos extra

K_ox=0.02;
K_f=0.2;
r=0.2;
