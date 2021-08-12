clear all; close all; clc

Trabajo_datos % Llamar a los datos del archivo Trabajo_datos
%% PRECAMARA
OF_st = 8;
% Calculo de temperatura de combustion
Q =  -(1 + 1/OF_st)*Entalp_h2o;

% Seleccionar el tipo de mezcla
syms OF_pc_sym
m_pc = 1;       % Seleccionar el tipo de mezcla en la precámara: 1 para rica (exceso de H2), 0 para pobre (exceso de O2)
if m_pc == 0    % Pobre --> Exceso de Oxidante
    Cp_pc = ( Cp_h2o + (OF_pc_sym/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OF_pc_sym/OF_st/2 - 1/2) );
    %Cp_pc = ( Cp_h2o_kg + (OF_pc_sym/OF_st/2 - 1/2)*Cp_o2_kg )/( 1 + (OF_pc_sym/OF_st/2 - 1/2) );
    Mp_pc = ( M_h2o + (OF_pc_sym/OF_st/2 - 1/2)*M_o2)/( 1 + (OF_pc_sym/OF_st/2 - 1/2) );
    % cp promedio en kg
    Cp_pc = Cp_pc/(Mp_pc*1e-3); %J/kgK
    fun = (T_tin-T_ref) == (Q/Cp_pc)*(OF_st/(1+OF_pc_sym));
else            % Rica --> Exceso de Fuel
    Cp_pc = ( Cp_h2o + (OF_st/OF_pc_sym - 1)*Cp_h2 )/( 1 + (OF_st/OF_pc_sym - 1) );
    Mp_pc = ( M_h2o + (OF_st./OF_pc_sym - 1)*M_h2 )/( OF_st/OF_pc_sym );
    % cp promedio en kg
    Cp_pc = Cp_pc/(Mp_pc*1e-3); %J/kgK
    fun = (T_tin == Q/Cp_pc*OF_pc_sym/(1.+OF_pc_sym) + T_ref);
end

OF_pc = max(double(vpasolve(fun)));  % O/F PC

% Repito para no simbolico
if m_pc == 0    % Pobre --> Exceso de Oxidante
        Cp_pc = ( Cp_h2o + (OF_pc/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OF_pc/OF_st/2 - 1/2) );
    Mp_pc = ( M_h2o + (OF_pc/OF_st/2 - 1/2)*M_o2)/( 1 + (OF_pc/OF_st/2 - 1/2) );
    % cp promedio en kg
    Cp_pc = Cp_pc/(Mp_pc*1e-3);
else            % Rica --> Exceso de Fuel
    Cp_pc = ( Cp_h2o + (OF_st/OF_pc - 1)*Cp_h2 )/( 1 + (OF_st/OF_pc - 1) );
    Mp_pc = ( M_h2o + (OF_st./OF_pc - 1)*M_h2 )/( OF_st/OF_pc );
    % cp promedio en kg
    Cp_pc = Cp_pc/(Mp_pc*1e-3);
end


% Valores necesarios de los productos:
R_pc = (R_un / (Mp_pc*1e-3)); %[J/kgK]
gamma_pc = Cp_pc / (Cp_pc - R_pc);
gamma_g_pc = sqrt(gamma_pc) * ((2 / (gamma_pc + 1))^((gamma_pc+1)/(2*(gamma_pc-1))));  % Gamma(gamma) de la pc

%% CAMARA COMBUSTION
OF = linspace(0.5, 20, 15);
Pc = linspace(1e5, 60e6, 15);
mt_m_v = zeros(length(OF), length(Pc));
Isp_0 = zeros(length(OF), length(Pc));
Isp_0_aux = zeros(length(OF), length(Pc));
Isp = zeros(length(OF), length(Pc));
T_c_save = zeros(length(OF), length(Pc));
OF_cc = zeros(length(OF), length(Pc));
mt_m = zeros(length(OF), length(Pc));

for i = 1:length(OF)
%  for i = 1:1
    for j = 1:length(Pc)
% Gastos másicos en las bombas 
    m_bf_m = 1 / (1 + OF(i)); % gasto másico bomba fuel
    m_bo_m = OF(i)/ (OF(i) + 1); % gasto másico bomba oxidante
    
% % Deltas P en las bombas
    Delta_P_bo = (Pc(j) / Pi_iny_cc_ox) - P_dox;
    Delta_P_bf = (Pc(j) / (Pi_iny_cc_red * Pi_refr)) - P_df;
    
% % Escogemos el gasto másico en la turbina / gasto másico total (mt_m) como
% % incógnita y resolvemos
    %syms mt_m   
    %S = vpasolve(((m_bo_m) * ( Delta_P_bo / (rend_bo * rho_ox) + m_bf_m * (Delta_P_bf / (rend_bf * rho_f))) == rend_mec * rend_t * mt_m * Cp_pc * T_tin * ( 1 - (1 / Pi_t) ^ ( ( gamma_pc - 1) / gamma_pc ))), mt_m);
   
    tau_t = rend_t*Cp_pc* T_tin*(1-(1/Pi_t)^((gamma_pc-1)/gamma_pc));
    tau_bo = Delta_P_bo / (rend_bo * rho_ox);
    tau_bf = Delta_P_bf / (rend_bf * rho_f);
    mt_m(i,j) = (m_bo_m*tau_bo+m_bf_m*tau_bf)/(rend_mec*tau_t);

% Gastos másicos
    m_o_pc = mt_m(i,j)* (OF_pc / (OF_pc + 1 )); %  gasto masico oxidante pre- camara / m total
    m_f_pc = mt_m(i,j) / (OF_pc + 1);            %  gasto masico fuel pre- camara / m total
    
    %m_o_c(i,j) = m_bo - m_o_pc;         %  gasto masico oxidante camara / m total
    m_o_c = m_bo_m - m_o_pc;         %  gasto masico oxidante camara / m total
    %m_f_c(i,j) = m_bf - m_f_pc;         %  gasto masico fuel camara / m total
    m_f_c = m_bf_m - m_f_pc;        %  gasto masico oxidante camara / m total
    
  % Reacción Química en la CC
    OF_cc(i,j) = (m_o_c/m_f_c);
     x_cc = OF_cc(i,j)*M_h2/M_o2;
     %A = m_f_c / 2;
     %B = m_o_c / 32;
     A = 1;
     B = x_cc;
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
   
    M_m = ( E * M_h2o + (C * M_h2) + (D * M_o2)) / (C+D+E); %[g/mol] Masa molecular media de los productos de la cc
    R_c = R_un*1000/ M_m; %[J/kgK]
    Cp_medio = ( E * Cp_h2o + (C * Cp_h2) + (D * Cp_o2)) / (C+D+E);  % [J/molK] Cp media de los productos de la cc 
    %Cp_medio_kg = (E*Cp_h2o_kg + C*Cp_h2_kg + D*Cp_o2_kg)/(E + C + D); % J/(g*K)
    Cp_medio_kg = Cp_medio/(M_m*10^(-3));
    
    gamma_c = 1 / (1 - (R_c / Cp_medio_kg)); % Gamma de la cc
    gamma_g = sqrt(gamma_c) * ((2 / (gamma_c + 1))^((gamma_c+1)/(2*(gamma_c-1))));  % Gamma(gamma) de la cc
    
    Entalp_h2o_mol = Entalp_h2o*(M_h2o);
        % PUEDE SER QUE HAYA QUE REDEFINIR LA
        % Q????????????????????????????????
    T_c = Q/Cp_medio_kg*OF_cc(i,j)/(1.+OF_cc(i,j)) + T_ref;
    c_star = sqrt(R_c * T_c)/ gamma_g;  % m/s
    
    func_p =  (2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc(j))^((gamma_c-1)/gamma_c));

   epsilon = gamma_g / (((P_s/Pc(j)) ^ (1 / gamma_c)) * sqrt(func_p));
   %epsilon = gamma_g / (((P_s/Pc(j)) ^ (1 / gamma_c)) * sqrt( (2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc(j))^((gamma_c-1)/gamma_c))));
   C_E = gamma_g* sqrt(func_p) + epsilon*(((P_s/Pc(j))) - (P_amb / Pc(j)));
   % Filtrar ISP (para eliminar datos con valores imaginarios o negativos
   % que no serían posibles)
   Isp_0(i,j) = C_E * c_star;
  
    %% TURBINA / TOBERA AUXILIAR
    
    % Calculo de presiones 
    P_e_t = (Pc(j) * Pi_iny_pc_ox * Pi_pc) / ( Pi_iny_cc_ox ); % P entrada turbina
    P_s_t = P_e_t / Pi_t; % P salida turbina   
    
    % Calculo parametros turbina
    T_s_t = T_tin * (1 - (rend_t*(1- Pi_t ^ ((gamma_pc-1)  / -gamma_pc))));  % T salida turbina
    
    func_RT = R_pc * T_s_t;
%     if func_RT < 0
%         func_RT = 0;
%     else
%     end
    c_star_aux = sqrt(func_RT)/ gamma_g_pc;
    ps_aux = P_s_t / ( 1 + (gamma_pc-1)/2)^(gamma_pc/(gamma_pc-1)); %  presion a la salida de le tobera auxiliar
    func_p_aux = (2*gamma_pc/(gamma_pc-1)) * (1 - (ps_aux/P_s_t)^((gamma_pc-1)/gamma_pc));
%     if func_p_aux < 0
%         func_p_aux = 0;
%     else
%     end
    %epsilon_aux = gamma_g_pc / (((ps_aux/P_s _t) ^ (1 / gamma_pc)) * sqrt( (2*gamma_pc/(gamma_pc-1)) * (1 - ((ps_aux/P_s_t))^((gamma_pc-1)/gamma_pc))));
    epsilon_aux = gamma_g_pc / (((ps_aux/P_s_t) ^ (1 / gamma_pc)) * sqrt(func_p_aux));
    %C_E_aux = gamma_g_pc * sqrt( (2*gamma_pc/(gamma_pc-1)) * (1 - (ps_aux/P_s_t)^((gamma_pc-1)/gamma_pc))) + epsilon_aux*(((ps_aux/P_s_t)) - (P_amb / P_s_t));
    C_E_aux = gamma_g_pc * sqrt(func_p_aux) + epsilon_aux*(((ps_aux/P_s_t)) - (P_amb / P_s_t));

    % FILTER ISP AUX

    Isp_0_aux(i,j) = C_E_aux * c_star_aux;
    
    %% ISP Total
 
    Isp = Isp_0 * (1 - mt_m(i,j)) + Isp_0_aux * mt_m(i,j); % Isp total
    %Isp(i,j) = Isp_0(i,j) * (1 - mt_m_v(i,j)) + Isp_0_aux(i,j) * mt_m_v(i,j);
     if ((mt_m(i,j)) > 1)  
        Isp(i,j) = 0; 
    else
        Isp(i,j) = Isp(i,j);
     end
    
    end
    
end
[i, j] = find(ismember(Isp, max(Isp(:))))

vx = Pc;              % Pc
vy = OF;               % O/F general
vo = OF_cc;


figure()
plot(OF_cc,Isp, 'DisplayName','OF_cc')
title('Isp vs O/F (CC)')
xlabel('O/F (CC)')
ylabel('Isp [s]')
axis tight
% 
figure()
plot(OF,Isp, 'DisplayName','OF')
title('Isp vs O/F')
xlabel('O/F')
ylabel('Isp [s]')
axis tight

figure()
surf(vx,vy,Isp)
title('Isp')
xlabel('Pc [Pa]')
ylabel('O/F')
zlabel('Isp [m/s]')
colorbar
axis tight
