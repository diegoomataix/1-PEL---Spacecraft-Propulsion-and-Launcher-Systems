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
    Cp_pc = ( Cp_h20 + (OF_pc_sym/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OF_pc_sym/OF_st/2 - 1/2) );
    Mp_pc = ( M_h20 + (OF_pc_sym/OF_st/2 - 1/2)*M_o2)/( 1 + (OF_pc_sym/OF_st/2 - 1/2) );
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
        Cp_pc = ( Cp_h20 + (OF_pc/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OF_pc/OF_st/2 - 1/2) );
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

OF = 8.3574;
Pc = 40e6;

% Gastos másicos en las bombas 
    m_bf_m = 1 / (1 + OF); % gasto másico bomba fuel
    m_bo_m = OF/ (OF + 1); % gasto másico bomba oxidante
    
% % Deltas P en las bombas
    Delta_P_bo = (Pc / Pi_iny_cc_ox) - P_dox;
    Delta_P_bf = (Pc / (Pi_iny_cc_red * Pi_refr)) - P_df;
    
% % Buscamos taus
    tau_t = rend_t*Cp_pc* T_tin*(1-(1/Pi_t)^((gamma_pc-1)/gamma_pc));
    tau_bo = Delta_P_bo / (rend_bo * rho_ox);
    tau_bf = Delta_P_bf / (rend_bf * rho_f);
    mt_m = (m_bo_m*tau_bo+m_bf_m*tau_bf)/(rend_mec*tau_t);

% Gastos másicos
    m_o_pc = mt_m* (OF_pc / (OF_pc + 1 )); %  gasto masico oxidante pre- camara / m total
    m_f_pc = mt_m / (OF_pc + 1);            %  gasto masico fuel pre- camara / m total
    
    m_o_c = m_bo_m - m_o_pc;         %  gasto masico oxidante camara / m total
    m_f_c = m_bf_m - m_f_pc;        %  gasto masico oxidante camara / m total
%    
%     
%     Q = -(1 + 1/OF_st)*Entalp_h2o;
%     OFr = 0.7477;
% %     cpp = ( Cp_h2o + (OF_st/OFp - 1)*Cp_h2 )/( 1 + (OF_st/OFp - 1) );
%     Cpr = ( Cp_h2o + (OFr/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OFr/OF_st/2 - 1/2) );
% 
%     % Masa molar promedio
% %     Mpp = ( Mh2o + (OFst./OFp - 1)*MH2 )./( OFst./OFp );
%     Mpr = ( M_h2o + (OFr/OF_st/2 - 1/2)*M_o2)./( 1 + (OFr/OF_st/2 - 1/2) );
% 
%     % cp promedio en kg
% %     cpp = cpp./(Mpp*1e-3);  
%     Cpr = Cpr/(Mpr*1e-3);
% 
%     % Temperatura de combustion
% %     Tcp = Q./cpp.*OFp./(1.+OFp) + Tref;
%     Tcr = Q/Cpr*OF_st/(1+OFr) + T_ref;
% 
%     % gamma promedio en kg
% %     gamma_p = cpp./(cpp - R./Mpp);
%     gamma_r = Cpr./(Cpr - R_un*1000/Mpr);
% 
%     % Propiedades de los productos
%     OF = OFr';
%     Cp_prod = [ Cpr];
%     M_prod = [ Mpr];
%     Tc_prod = [ Tcr];
%     gamma_c = [ gamma_r];
%     gamma_g = sqrt(gamma_c) * ((2 / (gamma_c + 1))^((gamma_c+1)/(2*(gamma_c-1))));
% 
%     % c* en funcion OF
%     c_star = sqrt(R_un/M_prod.*Tc_prod)./gamma_g;

    
  % Reacción Química en la CC
    OF_cc = (m_o_c/m_f_c);

    x_cc = OF_cc*M_h2/M_o2;
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
    Cp_medio_kg = Cp_medio/(M_m*10^(-3));
    
    gamma_c = 1 / (1 - (R_c / Cp_medio_kg)); % Gamma de la cc
    gamma_g = sqrt(gamma_c) * ((2 / (gamma_c + 1))^((gamma_c+1)/(2*(gamma_c-1))));  % Gamma(gamma) de la cc
    
    Entalp_h2o_mol = Entalp_h2o*(M_h2o);

    T_c = Q/Cp_medio_kg*OF_cc/(1.+OF_cc) + T_ref;
    c_star = sqrt(R_c * T_c)/ gamma_g;  % m/s
    
    func_p =  (2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc)^((gamma_c-1)/gamma_c));

    epsilon = gamma_c / (((P_s/Pc) ^ (1 / gamma_c)) * sqrt(func_p));
    C_E = gamma_c* sqrt(func_p) + epsilon*(((P_s/Pc)) - (P_amb / Pc));

    Isp_0 = C_E * c_star;
  
    %% TURBINA / TOBERA AUXILIAR
    
    % Calculo de presiones 
    P_e_t = (Pc * Pi_iny_pc_ox * Pi_pc) / ( Pi_iny_cc_ox ); % P entrada turbina
    P_s_t = P_e_t / Pi_t; % P salida turbina   
    
    % Calculo parametros turbina
    T_s_t = T_tin * (1 - (rend_t*(1- Pi_t ^ ((gamma_pc-1)  / -gamma_pc))));  % T salida turbina
    
    func_RT = R_pc * T_s_t;

    c_star_aux = sqrt(func_RT)/ gamma_g_pc;
    ps_aux = P_s_t / ( 1 + (gamma_pc-1)/2)^(gamma_pc/(gamma_pc-1)); %  presion a la salida de le tobera auxiliar
    func_p_aux = (2*gamma_pc/(gamma_pc-1)) * (1 - (ps_aux/P_s_t)^((gamma_pc-1)/gamma_pc));

    epsilon_aux = gamma_g_pc / (((ps_aux/P_s_t) ^ (1 / gamma_pc)) * sqrt(func_p_aux));
    C_E_aux = gamma_g_pc * sqrt(func_p_aux) + epsilon_aux*(((ps_aux/P_s_t)) - (P_amb / P_s_t));

    % FILTER ISP AUX

    Isp_0_aux = C_E_aux * c_star_aux;
    
    %% ISP Total
 
    Isp = Isp_0 * (1 - mt_m) + Isp_0_aux * mt_m % Isp total
    
    Isp_0
    
    Isp_0_aux
