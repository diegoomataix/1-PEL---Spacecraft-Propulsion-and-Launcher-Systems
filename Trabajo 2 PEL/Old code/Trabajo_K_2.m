clear all; close all; clc

Trabajo_datos % Llamar a los datos del archivo Trabajo_datos

%% PRECAMARA / TURBINA

% Seleccionar el tipo de mezcla
syms x_pc_sym
m_pc = 1;       % Seleccionar el tipo de mezcla en la precámara: 0 para rica (exceso de H2), 1 para pobre (exceso de O2)
if m_pc == 1    % Pobre --> Exceso de Oxidante
    E_pc = 1;
    C_pc = 0;
    D_pc = x_pc_sym - 1/2;
    
    M_prod_pc_sym = (E_pc*M_h20 + C_pc*M_h2 + D_pc*M_o2)/(E_pc + C_pc + D_pc); % kg/mol
    Cp_prod_pc_sym = (E_pc*Cp_h20 + C_pc*Cp_h2 + D_pc*Cp_o2)/((E_pc + C_pc + D_pc)); % J/(mol*K)
    Coef = (1+x_pc_sym)/8;

    Tc_pc = T_tin;

    Q = Coef * Cp_prod_pc_sym * (Tc_pc - T_ref) == -( E_pc / (E_pc + D_pc + C_pc)) * Entalp_h2o;  % J/mol
    x_pc_vect = double(vpasolve(Q));
    
else            % Rica --> Exceso de Fuel
    E_pc = 2*x_pc_sym;
    C_pc = 1-2*x_pc_sym;
    D_pc = 0;
    
    Cp_prod_pc_sym = (E_pc*Cp_h20 + C_pc*Cp_h2 + D_pc*Cp_o2)/((E_pc + C_pc + D_pc)); % J/(mol*K)
    Coef = (1+x_pc_sym)/x_pc_sym;
    
    Tc_pc = T_tin;

    Q = Coef * Cp_prod_pc_sym * (Tc_pc - T_ref) == -( E_pc / (E_pc + D_pc + C_pc)) * Entalp_h2o;  % J/mol
    x_pc_vect = double(vpasolve(Q));
end

x_pc = max(x_pc_vect(:));

% Para evitar que se guarden como sym usamos este otro if else
if m_pc == 1    % Pobre --> Exceso de Oxidante
    E_pc = 1;
    C_pc = 0;
    D_pc = x_pc - 1/2;
    
else            % Rica --> Exceso de Fuel
    E_pc = 2*x_pc;
    C_pc = 1-2*x_pc;
    D_pc = 0;
end

M_prod_pc = (E_pc*M_h20 + C_pc*M_h2 + D_pc*M_o2)/(E_pc + C_pc + D_pc); % kg/mol
Cp_prod_pc = (E_pc*Cp_h20 + C_pc*Cp_h2 + D_pc*Cp_o2)/((E_pc + C_pc + D_pc)*M_prod_pc); % J/(g*K)
Cp_prod_pc = Cp_prod_pc*1e3; % J/(kg*K)

OF_pc = x_pc * M_o2 / M_h2;
OF_pc_optimo = vpa(subs(OF_pc));% Usamos la funcion subs para sustituir en OF_pc el valor de x

% Valores necesarios de los productos:
R_pc = R_un *1000/ M_prod_pc; %[J/kgK]
gamma_pc = Cp_prod_pc / (Cp_prod_pc - R_pc);
gamma_g_pc = sqrt(gamma_pc) * ((2 / (gamma_pc + 1))^((gamma_pc+1)/(2*(gamma_pc-1))));  % Gamma(gamma) de la pc

%% EC ACOPLAMIENTO

% Definimos dos parámetros:
%   x = O/F (lo que se optimiza)
%   mt_m = m turbina / m total
% Les asignamos un rango de valores que se encuentre dentro de lo posible
% para un motor generador de gas. i.e.: max 100 bar Pc

% SET MESH

%x =  linspace(1, 20, 75);   % para ver o/f>=1
x =  linspace(6, 20, 10); % para ver o/f<=1
Pc = linspace(1e5, 50e6, 10);

% INITIALISE MATRICES

mt_m_v = zeros(length(x), length(Pc));
% m_o_c = zeros(length(x), length(Pc));
% m_f_c = zeros(length(x), length(Pc));
%C_E = zeros(length(x), length(Pc));
Isp_0 = zeros(length(x), length(Pc));
Isp_0_aux = zeros(length(x), length(Pc));
Isp = zeros(length(x), length(Pc));

% Usamos un loop para repetir la ecuación de acoplamiento, la reacción
% quimica en la CC y calcular el Isp para todo el rango de valores de Pc y
% O/F escogido

for i = 1:length(x)
    for j = 1:length(Pc)
% Gastos másicos en las bombas       
    m_bf = 1 / (1 + x(i)); % gasto másico bomba fuel
    m_bo = x(i) / (x(i) + 1); % gasto másico bomba oxidante
% Deltas P en las bombas
    Delta_P_bo = (Pc(j) / Pi_iny_cc_ox) - P_dox;
    Delta_P_bf = (Pc(j) / (Pi_iny_cc_red * Pi_refr)) - P_df;
% Escogemos el gasto másico en la turbina / gasto másico total (mt_m) como
% incógnita y resolvemos
    syms mt_m   
    S = vpasolve(((m_bo) * ( Delta_P_bo / (rend_bo * rho_ox) + m_bf * (Delta_P_bf / (rend_bf * rho_f))) == rend_mec * rend_t * mt_m * Cp_prod_pc * T_tin * ( 1 - (1 / Pi_t) ^ ( ( gamma_pc - 1) / gamma_pc ))), mt_m);

    mt_m_v(i, j) = S;

% Gastos másicos
    m_o_pc = S * (OF_pc / (OF_pc +1 )); %  gasto masico oxidante pre- camara / m total
    m_f_pc = S / (OF_pc + 1);            %  gasto masico fuel pre- camara / m total
    
    %m_o_c(i,j) = m_bo - m_o_pc;         %  gasto masico oxidante camara / m total
    m_o_c = m_bo - m_o_pc;         %  gasto masico oxidante camara / m total
    %m_f_c(i,j) = m_bf - m_f_pc;         %  gasto masico fuel camara / m total
    m_f_c = m_bf - m_f_pc;        %  gasto masico oxidante camara / m total
    
  % Reacción Química en la CC
    A = m_f_c / 2;
    B = m_o_c / 32;
    E = min(A, 2*B);
    
    if A > E
        C = A - E;
    else
        C = 0;
    end 
        
    if 2*B > E
            D = B - (1/2)* E;
    else
        D = 0;
    end
   
    M_m = ( E * M_h20 + (C * M_h2) + (D * M_o2)) / (C+D+E);    % Masa molecular media de los productos de la cc    
     
    R_c = R_un *1000/ M_m; %[J/kgK]
    
    Cp_medio = ( E * Cp_h20 + (C * Cp_h2) + (D * Cp_o2)) / (C+D+E);  % [J/molK] Cp media de los productos de la cc 
    Cp_medio_kg = Cp_medio*1000/M_m; %[J/kgK]
    
    gamma_c = 1 / (1 - (R_c / Cp_medio_kg)); % Gamma de la cc
    
    gamma_g = sqrt(gamma_c) * ((2 / (gamma_c + 1))^((gamma_c+1)/(2*(gamma_c-1))));  % Gamma(gamma) de la cc
    
    calor_h20=-Entalp_h2o*E/(C+D+E); % [J/mol]
    
%     T_c = calor_h20 * min(8, (m_o_c/m_f_c)) / ((1 + (m_o_c/m_f_c)) * Cp_medio) + T_ref;
    
    T_c = calor_h20 * min(8, (m_o_c/m_f_c)) / ((1 + (m_o_c/m_f_c)) * Cp_medio) + T_ref;
    
    c_star = sqrt(R_c * T_c)/ gamma_g; 
    epsilon = gamma_g / (((P_s/Pc(j)) ^ (1 / gamma_c)) * sqrt( (2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc(j))^((gamma_c-1)/gamma_c))));
    
    C_E = gamma_g* sqrt( (2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc(j)))^((gamma_c-1)/gamma_c) + epsilon*(((P_s/Pc(j))) - (P_amb / Pc(j))));
   % Filtrar ISP (para eliminar datos con valores imaginarios o negativos
   % que no serían posibles)
     
    if (class(c_star) == "complex") && (class(C_E) == "complex") 
        Isp_0(i,j) = 0; 
    else
        Isp_0(i,j) = C_E * c_star;
        
        if Isp_0(i,j) <= 0
            Isp_0(i,j) = 0;
        else
            Isp_0(i,j) = Isp_0(i,j);

        end
    end
    
    if mt_m_v(i,j)>1
    Isp_0(i,j) = 0;
    end
    if mt_m_v(i,j)<0
     Isp_0(i,j)=0;
    end
        
    %% TURBINA / TOBERA AUXILIAR
    
    % Calculo de presiones 
    P_e_t = (Pc(j) * Pi_iny_pc_ox * Pi_pc) / ( Pi_iny_cc_ox ); % P entrada turbina
    P_s_t = P_e_t / Pi_t; % P salida turbina   
    
    % Calculo parametros turbina
    T_s_t = T_tin * (1 - (rend_t*(1- Pi_t ^ ((gamma_pc-1)  / -gamma_pc))));  % T salida turbina
    
    c_star_aux = sqrt(R_pc * T_s_t)/ gamma_g_pc ;
    
    ps_aux = P_s_t / ( 1 + (gamma_pc-1)/2)^(gamma_pc/(gamma_pc-1)); %  presion a la salida de le tobera auxiliar
    
    %epsilon_aux = gamma_g_pc / (((ps_aux/P_s_t) ^ (1 / gamma_pc)) * sqrt( (2*gamma_pc/(gamma_pc-1)) * (1 - (ps_aux/P_s_t))^((gamma_pc-1)/gamma_pc)));
    epsilon_aux = 1;
    C_E_aux = gamma_g_pc * sqrt( (2*gamma_pc/(gamma_pc-1)) * (1 - (ps_aux/P_s_t))^((gamma_pc-1)/gamma_pc) + epsilon_aux*(((ps_aux/P_s_t)) - (P_amb / P_s_t)));
    %C_E_aux = gamma_g_pc * sqrt( (2*gamma_pc/(gamma_pc-1)) * (1 - (ps_aux/P_s_t))^((gamma_pc-1)/gamma_pc) + epsilon_aux*(((ps_aux/Pc(j))) - (P_amb / Pc(j))));
    
    % FILTER ISP AUX
    
    if (class(c_star_aux) == "complex") && (class(C_E_aux) == "complex") 
        Isp_0_aux(i,j) = 0; 
    else
         Isp_0_aux(i,j) = C_E_aux * c_star_aux;
        
        if Isp_0_aux(i,j) <= 0
            Isp_0_aux(i,j) = 0;
        else
            Isp_0_aux(i,j) = Isp_0_aux(i,j);

        end
    end
    if mt_m_v(i,j)>1
        Isp_0_aux(i,j) = 0;
    end
    if mt_m_v(i,j)<0
     Isp_0_aux(i,j)=0;
    end
    
    %% ISP Total
 
    Isp(i,j) = Isp_0(i,j) * (1 - mt_m_v(i,j)) + Isp_0_aux(i,j) * mt_m_v(i,j);    % Isp total
    
    if (class(mt_m_v(i,j)) == "complex")  
        Isp(i,j) = 0; 
    else
        Isp(i,j) = Isp(i,j);
        
    if mt_m_v(i,j)>=1
    Isp(i,j)=0;
    elseif mt_m_v(i,j)<0
    Isp(i,j)=0;
    end
    end
    end
    
end
[i, j] = find(ismember(Isp, max(Isp(:))))   % Find max Isp coordinates

%% 3D PLOTS
vx = Pc;           % Pc
vy = x;            % O/F

figure(1)
surf(vx,vy,Isp)
title('Isp')
xlabel('Pc [Pa]')
ylabel('O/F')
zlabel('Isp [m/s]')
colorbar
axis tight

figure(2)
surf(vx,vy,Isp_0)
title('Isp_0')
xlabel('Pc [Pa]')
ylabel('O/F')
zlabel('Isp [m/s]')
colorbar
axis tight

figure(3)
surf(vx,vy,Isp_0_aux)
title('Isp_0_(aux)')
xlabel('Pc [Pa]')
ylabel('O/F')
zlabel('Isp [m/s]')
colorbar
axis tight