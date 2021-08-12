clear all; close all; clc

Trabajo_datos % Llamar a los datos del archivo Trabajo_datos
%% PRECAMARA
OF_st = 8;
% Calculo de temperatura de combustion
Q =  -(1 + 1/OF_st)*Entalp_h2o;

% Seleccionar el tipo de mezcla
syms OF_pc_sym
m_pc = 1;       % Seleccionar el tipo de mezcla en la precÃ¡mara: 0 para rica (exceso de H2), 1 para pobre (exceso de O2)
if m_pc == 0    % Pobre --> Exceso de Oxidante
    Cp_pc = ( Cp_h20 + (OF_pc_sym/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OF_pc_sym/OF_st/2 - 1/2) );
    Mp_pc = ( M_h20 + (OF_pc_sym/OF_st/2 - 1/2)*M_o2)/( 1 + (OF_pc_sym/OF_st/2 - 1/2) );
    % cp promedio en kg
    Cp_pc = Cp_pc/(Mp_pc*1e-3); %J/kgK
    fun = (T_tin-T_ref) == (Q/Cp_pc)*(OF_st/(1+OF_pc_sym));
else            % Rica --> Exceso de Fuel
    Cp_pc = ( Cp_h20 + (OF_st/OF_pc_sym - 1)*Cp_h2 )/( 1 + (OF_st/OF_pc_sym - 1) );
    Mp_pc = ( M_h20 + (OF_st./OF_pc_sym - 1)*M_h2 )/( OF_st/OF_pc_sym );
    % cp promedio en kg
    Cp_pc = Cp_pc/(Mp_pc*1e-3); %J/kgK
    fun = (T_tin == Q/Cp_pc*OF_pc_sym/(1.+OF_pc_sym) + T_ref);
end

OF_pc = max(double(vpasolve(fun)));  % O/F PC

% Repito para no simbolico
if m_pc == 0    % Pobre --> Exceso de Oxidante
        Cp_pc = ( Cp_h20 + (OF_pc/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OF_pc/OF_st/2 - 1/2) );
    Mp_pc = ( M_h20 + (OF_pc/OF_st/2 - 1/2)*M_o2)/( 1 + (OF_pc/OF_st/2 - 1/2) );
    % cp promedio en kg
    Cp_pc = Cp_pc/(Mp_pc*1e-3);
else            % Rica --> Exceso de Fuel
    Cp_pc = ( Cp_h20 + (OF_st/OF_pc - 1)*Cp_h2 )/( 1 + (OF_st/OF_pc - 1) );
    Mp_pc = ( M_h20 + (OF_st./OF_pc - 1)*M_h2 )/( OF_st/OF_pc );
    % cp promedio en kg
    Cp_pc = Cp_pc/(Mp_pc*1e-3);
end

% Valores de los productos:
R_pc = (R_un / (Mp_pc*1e-3)); %[J/kgK]
gamma_pc = Cp_pc / (Cp_pc - R_pc);
gamma_g_pc = sqrt(gamma_pc) * ((2 / (gamma_pc + 1))^((gamma_pc+1)/(2*(gamma_pc-1))));  % Gamma(gamma) de la pc

%% CAMARA COMBUSTION
OF_cc = linspace(0.5, 12, 50);
Pc = linspace(1e5, 50e6, 50);

 Isp_0 = zeros(length(OF_cc), length(Pc));
 Isp_0_aux = zeros(length(OF_cc), length(Pc));
 Isp = zeros(length(OF_cc), length(Pc));
 T_c = zeros(length(OF_cc), length(Pc));
 c_star = zeros(length(OF_cc), length(Pc));
 C_E = zeros(length(OF_cc), length(Pc));
 C_E_aux = zeros(length(OF_cc), length(Pc));
 K_m = zeros(length(OF_cc), length(Pc));

for i = 1:length(OF_cc)
    for j = 1:length(Pc)
        
% Gastos másicos en las bombas:
%    m_bf_m = 1 / (1 + OF(i)); % gasto másico bomba fuel
%    m_bo_m = OF(i)/ (OF(i) + 1); % gasto másico bomba oxidante
    
% % Deltas P en las bombas:
    Delta_P_bo = (Pc(j) / Pi_iny_cc_ox) - P_dox;
    Delta_P_bf = (Pc(j) / (Pi_iny_cc_red * Pi_refr)) - P_df;
    
% Definismos los taus de las bombas y la turbina:
    tau_t = rend_t*Cp_pc* T_tin*(1-(1/Pi_t)^((gamma_pc-1)/gamma_pc));
    tau_bo = Delta_P_bo / (rend_bo * rho_ox);
    tau_bf = Delta_P_bf / (rend_bf * rho_f);
   
% Gastos másicos:
    if OF_cc(i) > 8 %Mezcla Pobre
        Cp_cc = ( Cp_h20 + (OF_cc(i)/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OF_cc(i)/OF_st/2 - 1/2) );
        Mp_cc = ( M_h20 + (OF_cc(i)/OF_st/2 - 1/2)*M_o2)/( 1 + (OF_cc(i)/OF_st/2 - 1/2) );
        % cp promedio en kg
        Cp_cc = Cp_cc/(Mp_cc*1e-3); %J/kgK
        T_c(i,j) = (Q/Cp_cc)*(OF_st/(1+OF_cc(i))) + T_ref;
    else 
        Cp_cc = ( Cp_h20 + (OF_st/OF_cc(i) - 1)*Cp_h2 )/( 1 + (OF_st/OF_cc(i) - 1) );
        Mp_cc = ( M_h20 + (OF_st/OF_cc(i) - 1)*M_h2 )/( OF_st/OF_cc(i) );
        % cp promedio en kg
        Cp_cc = Cp_cc/(Mp_cc*1e-3); %J/kgK
        T_c(i,j) = Q/Cp_cc*OF_cc(i)/(1+OF_cc(i)) + T_ref;
    end
    
    R_c = R_un*1000/ Mp_cc; %[J/kgK]
    gamma_c = 1 / (1 - (R_c / Cp_cc)); % Gamma de la cc
    gamma_g = sqrt(gamma_c) * ((2 / (gamma_c + 1))^((gamma_c+1)/(2*(gamma_c-1))));  % Gamma(gamma) de la cc

    c_star(i,j) = sqrt(R_c * T_c(i,j))/ gamma_g;  % m/s
    
    func_p =  (2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc(j))^((gamma_c-1)/gamma_c));
    epsilon = gamma_g / (((P_s/Pc(j)) ^ (1 / gamma_c)) * sqrt(func_p));
    C_E(i,j) = gamma_g* sqrt(func_p) + epsilon*(((P_s/Pc(j))) - (P_amb / Pc(j)));
    Isp_0(i,j) = C_E(i,j) * c_star(i,j);
  
    %% TURBINA / TOBERA AUXILIAR
    
    % Calculo de presiones:
    P_e_t = (Pc(j) * Pi_iny_pc_ox * Pi_pc) / ( Pi_iny_cc_ox ); % P entrada turbina
    P_s_t = P_e_t / Pi_t; % P salida turbina   
    
    % Calculo parametros turbina:
    T_s_t = T_tin * (1 - (rend_t*(1- (1/Pi_t) ^ ((gamma_pc-1)  / gamma_pc))));  % T salida turbina 
    func_RT = R_pc * T_s_t;

    c_star_aux = sqrt(func_RT)/ gamma_g_pc;
    ps_aux = P_s_t / ( 1 + (gamma_pc-1)/2)^(gamma_pc/(gamma_pc-1)); %  presion a la salida de le tobera auxiliar
    func_p_aux = (2*gamma_pc/(gamma_pc-1)) * (1 - (ps_aux/P_s_t)^((gamma_pc-1)/gamma_pc));
    epsilon_aux = gamma_g_pc / (((ps_aux/P_s_t) ^ (1 / gamma_pc)) * sqrt(func_p_aux));
    C_E_aux(i,j) = gamma_g_pc * sqrt(func_p_aux) + epsilon_aux*(((ps_aux/P_s_t)) - (P_amb / P_s_t));

    % ISP AUX
    Isp_0_aux(i,j) = C_E_aux(i,j) * c_star_aux;
    
    %% Acoplamiento de la turbomaquinaria
    moC_moPC = ( (1 + 1/OF_pc)*rend_mec* tau_t - (tau_bo + tau_bf/OF_pc) )/...
               (tau_bo + tau_bf/OF_cc(i));
    mT_moPC = moC_moPC*(1 + 1/OF_cc(i)) + (1 + 1/OF_pc);

    %% ISP Total
 
    Isp(i,j) = Isp_0(i,j) * (1 -(1 + 1/OF_pc)/mT_moPC) + Isp_0_aux(i,j) * (1 + 1/OF_pc)/mT_moPC; % Isp total
    
%     %% Km
%     
%     K_m(i,j) = ((K_o * OF) + K_f) / ( 1 + OF);
%     
    end
    
end
[i, j] = find(ismember(Isp, max(Isp(:))))   

%% GASTO MÁSICO
m_t = E ./ Isp;

mo_PC = m_t ./ mT_moPC; 
mf_PC = mo_PC / OF_pc;

mo_C = mo_PC .* moC_moPC; 
mf_C = mo_C / OF_cc;

OF = (mo_C + mo_PC)./(mf_C + mf_PC);

%%  OPTIMIZACIÓN (APARTADO 2)

Isp_Km = Isp / (1+K_m);

%%  OPTIMIZACIÓN (APARTADO 2) - IMPULSO VOLUMÉTRICO

rho_m = (1 + OF) / ( (OF / rho_ox) + (1 / rho_f));

Isp_Vol = Isp(i,j) * rho_m;

%%  OPTIMIZACIÓN (APARTADO 2) - DELTA V

DV = Isp * log((1 + K_m)/(r + K_m));

%% 2D PLOTS

% Definir variables para pintar curvas

vx = round(Pc,-3);              % Pc
% vxr(:) = round(vx(:), 0);
vy = OF;                        % O/F general
vo = OF_cc;

%%%%%%%%%%%%%%%%%%%%%% % 2D PLOT / ISP vs OF (cc)
figure()    
plot(vo,Isp/g0, 'DisplayName','x')
hold on
title('Isp vs OF')
xlabel('O/F (Cámara de Combustión)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Isp [s]','FontSize',12,'FontWeight','bold','Color','k')
lgd = legend(strcat(strtrim(cellstr(num2str(single(vx(:)/1e5))))),'NumColumns',2);
grid on
title(lgd, 'Pc [bar]')
xl = xlim;
yl = ylim;
line([8, 8], yl); % Old way of doing xline().
axis tight
axis([0 12 2500/g0 4100/g0])
hold off

%%%%%%%%%%%%%%%%%%%%%% % 2D PLOT / ISP_0 vs OF (cc)
figure()   
plot(vo,Isp_0/g0, 'DisplayName','x')
hold on
title('Isp^0 vs OF')
xlabel('O/F (Cámara de Combustión)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Isp^0 [s]','FontSize',12,'FontWeight','bold','Color','k')
lgd = legend(strcat(strtrim(cellstr(num2str(single(vx(:)/1e5))))),'NumColumns',2);
grid on
title(lgd, 'Pc [bar]')
xl = xlim;
yl = ylim;
line([8, 8], yl); % Old way of doing xline().
axis tight
axis([0 12 2500/g0 4300/g0])
hold off

%%%%%%%%%%%%%%%%%%%%%%% 2D PLOT / ISP_0_aux vs OF (cc)
figure()    
plot(vo,Isp_0_aux/g0, 'DisplayName','x')
hold on
title('Isp^0_{aux} vs OF')
xlabel('O/F (Cámara de Combustión)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Isp^0_{aux} [s]','FontSize',12,'FontWeight','bold','Color','k')
lgd = legend(strcat(strtrim(cellstr(num2str(single(vx(:)/1e5))))),'NumColumns',2);
grid on
title(lgd, 'Pc [bar]')
% xl = xlim;
% yl = ylim;
% line([OF_pc, OF_pc], yl); % Old way of doing xline().
axis tight
axis([0 12 500/g0 2800/g0])
hold off

%%%%%%%%%%%%%%%%%%%%%%% 2D PLOT / Tc vs OF (cc)
figure()    
p1 = plot(vo,T_c, 'r', 'DisplayName','x');
p1(1).LineWidth = 2;
hold on
p2 = plot(vo,c_star, 'g','DisplayName','x');
p2(1).LineWidth = 2;
title('T_c & c^* vs OF')
xlabel('O/F (Cámara de Combustión)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('T_c [K], c^* [m/s]','FontSize',12,'FontWeight','bold','Color','k')
grid on
yl = ylim;
line([OF_st, OF_st], yl); % Old way of doing xline().
%legend('Tc','c^*', 'Estequiométrico');
axis tight
axis([0 12 500 5000])
hold off

%%%%%%%%%%%%%%%%%%%%%%% 2D PLOT / C_E vs OF (cc)
figure()    
% plot(vo,C_E, 'b','DisplayName','x')
plot(vo,C_E,'DisplayName','x')
hold on
%plot(vo,C_E_aux, 'm','DisplayName','x');
title('C_E vs OF')
xlabel('O/F (Cámara de Combustión)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('C_E','FontSize',12,'FontWeight','bold','Color','k')
%title('C_E & C_{E,aux} vs OF')
%xlabel('O/F (Cámara de Combustión)','FontSize',12,'FontWeight','bold','Color','k')
%ylabel('C_E','FontSize',12,'FontWeight','bold','Color','k')
grid on
xl = xlim;
yl = ylim;
line([OF_st, OF_st], yl); % Old way of doing xline().
axis tight
axis([0 12 .25 2])
hold off

%%%%%%%%%%%%%%%%%%%%%%% 2D PLOT / Gasto oxidante Cámara vs OF (cc)
figure()    
plot(vo,mo_C, 'DisplayName','x')
hold on
title('Gasto oxidante Cámara vs OF')
xlabel('O/F (Cámara de Combustión)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Gasto oxidante Cámara','FontSize',12,'FontWeight','bold','Color','k')
lgd = legend(strcat(strtrim(cellstr(num2str(single(vx(:)/1e5))))),'NumColumns',2);
grid on
yl = ylim;
line([OF_st, OF_st], yl); % Old way of doing xline().
title(lgd, 'Pc [bar]')
axis tight
axis([0 12 75 150])
hold off

%%%%%%%%%%%%%%%%%%%%%%% 2D PLOT / Gasto fuel y oxidante Precámara vs OF (cc)
figure()    
plot(vo,mf_PC, 'DisplayName','x')
hold on
plot(vo,mo_PC, 'DisplayName','x')
title('Gasto fuel Precámara vs OF')
xlabel('O/F (Cámara de Combustión)','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Gasto fuel Precámara','FontSize',12,'FontWeight','bold','Color','k')
lgd = legend(strcat(strtrim(cellstr(num2str(single(vx(:)/1e5))))),'NumColumns',2);
grid on
title(lgd, 'Pc [bar]')
axis tight
axis([0 12 0 75])
hold off

%% 3D PLOTS

%%%%%%%%%%%%%%%%%%%%%% 3D PLOT
figure()
surf(vo,(vx(:)/1e5),Isp/g0)
title('Isp')
xlabel('O/F (CC)')
ylabel('Pc [bar]')
zlabel('Isp [s]')
colorbar
axis tight
axis([0 12 2500/g0 4100/g0 0 410])

figure()
surf(vo,(vx(:)/1e5),T_c)
xlabel('O/F (CC)')
ylabel('Pc [bar]')
ylabel('O/F (CC)')
zlabel('T_c [K]')
colorbar
axis tight

figure()
surf(vo,(vx(:)/1e5),DV)
xlabel('O/F (CC)')
ylabel('Pc [bar]')
ylabel('O/F (CC)')
zlabel('DV [m/s]')
colorbar
axis tight

% 
% vect_c1 = linspace(max(max(Isp)),min(min(Isp)),20);
% vect_c1(20) = 3880.4;
% vect_c1(21) = 3920.4;
% vect_c1(22) = 3927.4;
% figure()
% hold on
% colormap(jet)
% grid on
% contour(vx, vo, Isp,[vect_c1(:)],'ShowText','on');


