clear all; close all; clc

%% DATOS

% Densidades [kg/m^3]
rho_ox = 1140; 
rho_f = 71;

% Densidades [Pa]
P_s = 40000; % Pa
P_amb = 101325; % Pa {SL}

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

% Presiones de los depósitos
P_dox = 6*10^5; % Pa
P_df = 3*10^5; %Pa

% Temperaturas
T_ref = 298.15;
T_pc = 900;
T_tin = 900;

% Cp [J/molK]
Cp_h20 = 54.32;
Cp_o2 = 39.07;
Cp_h2 = 36.11;

% Masas molares [kg/mol] y constantes
M_h20 = 18;
M_o2 = 32;
M_h2 = 2;

Entalp_h2o = -241818;
R_un = 8314.4598;

%% PRECAMARA / TURBINA

syms x_pc_sym
m_pc = 1;       % Seleccionar el tipo de mezcla en la precámara: 1 para pobre (exceso de O2), 0 para rica (exceso de H2)
if m_pc == 1    % Pobre --> Exceso de Oxidante
    E_pc = 1;
    C_pc = 0;
    D_pc = x_pc_sym - 1/2;
    Mezcla = 'pobre';
    
else            % Rica --> Exceso de Fuel
    E_pc = 2*x_pc_sym;
    C_pc = 1-2*x_pc_sym;
    D_pc = 0;
    Mezcla = 'rica';
end

M_prod_pc_sym = (E_pc*M_h20 + C_pc*M_h2 + D_pc*M_o2)/(E_pc + C_pc + D_pc); % kg/mol
Cp_prod_pc_sym = (E_pc*Cp_h20 + C_pc*Cp_h2 + D_pc*Cp_o2)/((E_pc + C_pc + D_pc)*M_prod_pc_sym); % J/(kg*K)

% syms Tc
Tc_pc = T_tin;

Q = M_prod_pc_sym * Cp_prod_pc_sym * (Tc_pc - T_ref) == -( E_pc / (E_pc + D_pc + C_pc)) * Entalp_h2o;  % J/mol
x_pc = double(vpasolve(Q));

if m_pc == 1    % Pobre --> Exceso de Oxidante
    E_pc = 1;
    C_pc = 0;
    D_pc = x_pc - 1/2;
    Mezcla = 'pobre';
    
else            % Rica --> Exceso de Fuel
    E_pc = 2*x_pc;
    C_pc = 1-2*x_pc;
    D_pc = 0;
    Mezcla = 'rica';
end

M_prod_pc = (E_pc*M_h20 + C_pc*M_h2 + D_pc*M_o2)/(E_pc + C_pc + D_pc); % kg/mol
Cp_prod_pc = (E_pc*Cp_h20 + C_pc*Cp_h2 + D_pc*Cp_o2)/((E_pc + C_pc + D_pc)*M_prod_pc); % J/(kg*K)

OF_pc = x_pc * M_o2 / M_h2;
OF_pc_optimo = vpa(subs(OF_pc));
% Subs para sustituir en OF el valor obtenido de x

% Valores necesarios de los productos:
R_pc = R_un / M_prod_pc;
gamma_pc = Cp_prod_pc / (Cp_prod_pc - R_pc);
gamma_g_pc = sqrt(gamma_pc) * ((2 / (gamma_pc + 1))^((gamma_pc+1)/(2*(gamma_pc-1))));  % Gamma(gamma) de la pc

%% EC ACOPLAMIENTO

% x = O/F (lo que se optimiza)
% mt_m = m turbina / m total

% SET MESH

x =  linspace(1, 50, 10);
Pc = linspace(1e5, 75e5, 10);

% INITIALISE MATRICES

mt_m_v = zeros(length(x), length(Pc));
% m_o_c = zeros(length(x), length(Pc));
% m_f_c = zeros(length(x), length(Pc));
%C_E = zeros(length(x), length(Pc));
Isp_0 = zeros(length(x), length(Pc));
Isp_0_aux = zeros(length(x), length(Pc));
Isp = zeros(length(x), length(Pc));

for i = 1:length(x)
    for j = 1:length(Pc)
    %for j = 1:1
    
    m_bf = 1 / (1 + x(i)); % gasto másico bomba fuel
    m_bo = x(i) / (x(i) + 1); % gasto másico bomba oxidante
    Delta_P_bo = (Pc(j) / Pi_iny_cc_ox) - P_dox;
    Delta_P_bf = (Pc(j) / (Pi_iny_cc_red * Pi_refr)) - P_df;

    syms mt_m

    % eqn = (m_bo) * ( Delta_P_bo / (rend_bo * rho_ox) + m_bf * (Delta_P_bf / (rend_bf * rho_f))) == rend_mec * rend_t * Cpp * T_tin * ( 1 - Pi_t ^ ( ( gamma - 1) / gamma ));

   % S = vpasolve(((m_bo) * ( Delta_P_bo / (rend_bo * rho_ox) + m_bf * (Delta_P_bf / (rend_bf * rho_f))) == rend_mec * rend_t * mt_m * Cpp * T_tin * ( 1 - Pi_t ^ ( ( gamma - 1) / gamma ))), mt_m)
   
    S = vpasolve(((m_bo) * ( Delta_P_bo / (rend_bo * rho_ox) + m_bf * (Delta_P_bf / (rend_bf * rho_f))) == rend_mec * rend_t * mt_m * Cp_prod_pc * T_tin * ( 1 - (1 / Pi_t) ^ ( ( gamma_pc - 1) / gamma_pc ))), mt_m);

    mt_m_v(i, j) = S;
    
    m_o_pc = S * (OF_pc / (OF_pc +1 )); %  gasto masico oxidante pre- camara / m total
    
    m_f_pc = S / (OF_pc + 1);            %  gasto masico fuel pre- camara / m total
    
    %m_o_c(i,j) = m_bo - m_o_pc;         %  gasto masico oxidante camara / m total
    m_o_c = m_bo - m_o_pc;         %  gasto masico oxidante camara / m total
    %m_f_c(i,j) = m_bf - m_f_pc;         %  gasto masico fuel camara / m total
    m_f_c = m_bf - m_f_pc;        %  gasto masico oxidante camara / m total
    
    
  % CHEMICAL REACTION
    
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
     
    R_c = R_un / M_m;
    
    Cp_medio = ( E * Cp_h20 + (C * Cp_h2) + (D * Cp_o2)) / (C+D+E);  % Cp media de los productos de la cc 
    
    gamma_c = 1 / (1 - (R_c / Cp_medio));           % Gamma de la cc
    
    gamma_g = sqrt(gamma_c) * ((2 / (gamma_c + 1))^((gamma_c+1)/(2*(gamma_c-1))));  % Gamma(gamma) de la cc
    
   % 
    calor_h20=-Entalp_h2o*E/(C+D+E);
    T_c = calor_h20 * min(8, (m_o_c/m_f_c)) / ((1 + (m_o_c/m_f_c)) * Cp_medio) + T_ref;
    
    c_star = sqrt(R_c * T_c)/ gamma_g; 
    
    epsilon = gamma_g / (((P_s/Pc(j)) ^ (1 / gamma_c)) * sqrt( (2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc(j))^((gamma_c-1)/gamma_c))));
    
    C_E = gamma_g * sqrt( (2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc(j)))^((gamma_c-1)/gamma_c)) + epsilon*((P_s/Pc(j)) - (P_amb / Pc(j)));
    
   % FILTER ISP
    
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
    
    %% TURBINA / TOBERA AUXILIAR
    
    P_e_t = (Pc(j) * Pi_iny_pc_ox * Pi_pc) / ( Pi_iny_cc_ox ); % P entrada turbina
    P_s_t = P_e_t / Pi_t; % P salida turbina   
    
    T_s_t = T_tin * (1 - (rend_t*(1- Pi_t ^ ((gamma_pc-1)  / -gamma_pc))));  % T salida turbina
    
    c_star_aux = sqrt(R_pc * T_s_t)/ gamma_g_pc ;
    
    ps_aux = P_s_t / ( 1 + (gamma_pc-1)/2)^(gamma_pc/(gamma_pc-1)); %  presion a la salida de le tobera auxiliar
    
    epsilon_aux = gamma_g_pc / (((ps_aux/P_s_t) ^ (1 / gamma_pc)) * sqrt( (2*gamma_pc/(gamma_pc-1)) * (1 - (ps_aux/P_s_t))^((gamma_pc-1)/gamma_pc)));
    
    C_E_aux = gamma_g_pc * sqrt( (2*gamma_pc/(gamma_pc-1)) * (1 - (ps_aux/P_s_t))^((gamma_pc-1)/gamma_pc) + epsilon_aux*(((ps_aux/P_s_t)) - (P_amb / P_s_t)));

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

% PLOT
figure(1)
plot(Pc, Isp_0)
xlabel("Pc")
ylabel("Isp_0")

figure(2)
plot(Pc, Isp_0_aux)
xlabel("Pc")
ylabel("Isp_0_aux")

figure(3)
plot(Pc, Isp)
xlabel("Pc")
ylabel("Isp")

figure(4)
plot(x, Isp_0)
xlabel("O/F")
ylabel("Isp_0")

figure(5)
plot(x, Isp_0_aux)
xlabel("O/F")
ylabel("Isp_0_aux")

figure(6)
plot(x, Isp)
xlabel("O/F")
ylabel("Isp")

[i, j] = find(ismember(Isp, max(Isp(:))));