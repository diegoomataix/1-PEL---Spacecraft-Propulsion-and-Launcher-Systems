
clear all, close all
%% Probl. 2

% Para la precámara, consideramos tanto mezclas ricas como pobres:
% H2 + x*O2 -> a*H2O + b*H2 + c*O2 + Q
% Mezcla rica: x<0.5 -> c=0; a=2x; b=1-2x
% Mezcla esteq: x = 0.5 -> a = 1; b = c = 0
% Mezcla pobre: x>0.5 -> a = 1; b = 0; c = x-1/2

Datos_probl2 % Llamamos a todos los valores datos del probl 2

%% Variación de parámetros
% Variación OF
OF_inicial = 0.1;
OF_final = 20;
n_puntos_OF = 20;
paso_OF = (OF_final - OF_inicial)/(n_puntos_OF - 1);

% Variación Pc [Pa]
Pc_inicial = 5*1e6; % 1MPa
Pc_final = 50*1e6; % 15MPa
n_puntos_Pc = 20;
paso_Pc = (Pc_final - Pc_inicial)/(n_puntos_Pc - 1);

%% Precámara
syms x_pc
m_pc = 0; % Tipo de mezcla en la precámara, 0 para pobre (exceso de O2), 1 para rica (exceso de H2)
if m_pc == 1 % Pobre -- Exceso de Oxidante
    a_pc = 1;
    b_pc = 0;
    c_pc = x_pc - 1/2;
    Mezcla = 'pobre';
else % Rica
    a_pc = 2*x_pc;
    b_pc = 1-2*x_pc;
    c_pc = 0;
    Mezcla = 'rica';
end

M_prod_pc = (a_pc*M_H2O + b_pc*M_H2 + c_pc*M_O2)/(a_pc + b_pc + c_pc); % kg/mol
Cp_prod_pc = (a_pc*Cp_H2O + b_pc*Cp_H2 + c_pc*Cp_O2)/((a_pc + b_pc + c_pc)*M_prod_pc); % J/(kg*K)
% Cp_prod_pc = Cp_prod_pc*1e3; % J/(kg*K)

% syms Tc
Tc_pc = T_max_t;
Q1_pc = (a_pc+b_pc+c_pc)*M_prod_pc*Cp_prod_pc*(Tc_pc - T_ref); % J
Q2_pc = -a_pc*Entalp_H2O; % J

OF_pc = x_pc*M_O2/M_H2;

fun_pc = Q1_pc == Q2_pc;
x_pc = solve(fun_pc);

OF_pc_optimo = vpa(subs(OF_pc));% vpa es solo para mostrarlo en formato decimal en vez de fraccion que muestra por defecto
% Subs para sustituir en OF el valor obtenido de x

% Valores necesarios de los productos:
R_pc = Ru/M_prod_pc;
gamma_pc = Cp_prod_pc/(Cp_prod_pc - R_pc);
Gamma_pc = sqrt(gamma_pc)*(2/(gamma_pc+1))^((gamma_pc+1)/(2*(gamma_pc-1)));

fprintf('Para mezcla %s, el O/F que maximiza la temperatura\n a la entrada de la turbina es: O/F = %.4f\n',Mezcla,OF_pc_optimo)

%% Sistema completo

syms Pc_cc


% Calculo de presiones

DeltaP_b_f = Pc_cc/(pi_ref*pi_iny_cc_f) - P_dep_f;
DeltaP_b_ox = Pc_cc/pi_iny_cc_ox - P_dep_ox;

P_iny_pc_f = P_dep_f + DeltaP_b_f;
P_iny_pc_ox = P_dep_ox + DeltaP_b_ox;

P_iny_cc_f = P_iny_pc_f*pi_ref;
P_iny_cc_ox = P_iny_pc_ox;

% Precamara
P_pc = P_iny_pc_ox*pi_iny_pc_ox;

P_entrada_turbina = P_pc*pi_pc;
P_salida_turbina = P_entrada_turbina*pi_t_max;

% Gastos másicos:
syms OF_cc

syms mf_pc mf_cc % de momento hasta q sepa calcularlas

mox_pc = mf_pc*OF_pc;
mox_cc = mf_cc*OF_cc;

mf = mf_cc + mf_pc;
mox = mox_cc + mox_pc;
mt = mf_pc + mox_pc;

% Todo lo anterior depende de Pc, mf_pc y mf_pc, necesito 2 ecuaciones más
% para hacerlo depender solo de Pc

% 1 - Acoplamiento turbina:
trab_turbina = rend_t*mt*Cp_prod_pc*T_max_t*(1 - (pi_t_max)^((gamma_pc-1)/gamma_pc));
trab_bomba_ox = mox*DeltaP_b_ox/(rend_b_ox*densidad_ox);
trab_bomba_f = mf*DeltaP_b_f/(rend_b_f*densidad_f);

fun1 = subs(rend_mec*trab_turbina == trab_bomba_ox + trab_bomba_f); % depende de Pc, mf_pc y mf_cc

syms a_cc b_cc c_cc
M_prod_cc = (a_cc*M_H2O + b_cc*M_H2 + c_cc*M_O2)/(a_cc + b_cc + c_cc); % kg/mol
Cp_prod_cc = (a_cc*Cp_H2O + b_cc*Cp_H2 + c_cc*Cp_O2)/((a_cc + b_cc + c_cc)*M_prod_cc); % J/(kg*K)
% Cp_prod_cc = Cp_prod_cc*1e3; % J/(kg*K)
    
%     Q1_cc = (a_cc+b_cc+c_cc)*M_prod_cc*Cp_prod_cc*(Tc_cc - T_ref); % J
%     Q2_cc = -a_cc*Entalp_H2O; % J
Tc_cc = T_ref - a_cc*Entalp_H2O/((a_cc+b_cc+c_cc)*M_prod_cc*Cp_prod_cc);

% Valores necesarios de los productos:
R_cc = Ru/M_prod_cc; % J/(molK/(g/mol)) = J/(K*g)
gamma_cc = Cp_prod_cc/(Cp_prod_cc - R_cc); % ambos están en J/g*K
Gamma_cc = sqrt(gamma_cc)*(2/(gamma_cc+1))^((gamma_cc+1)/(2*(gamma_cc-1)));

% E = mt*Isp_pc + (mox_cc + mf_cc)*Isp_cc
c_est_pc = sqrt(R_pc*1e3*Tc_pc)/Gamma_pc; 
c_est_cc = sqrt(R_cc*1e3*Tc_cc)/Gamma_cc;

Pc_pc = P_salida_turbina;
Ps_Pc_pc = P_s/Pc_pc;
Ps_Pc_cc = P_s/Pc_cc;
Pamb_Pc_pc = 101325/Pc_pc; % Pa/Pa
Pamb_Pc_cc = 101325/Pc_cc;
% eps_pc = Gamma_pc/(Ps_Pc_pc^(1/gamma_pc)*sqrt(2*gamma_pc/(gamma_pc-1)*(1-Ps_Pc_pc^((gamma_pc-1)/gamma_pc))));
eps_cc = Gamma_cc/(Ps_Pc_cc^(1/gamma_cc)*sqrt(2*gamma_cc/(gamma_cc-1)*(1-Ps_Pc_cc^((gamma_cc-1)/gamma_cc))));
C_E_pc = Gamma_pc*sqrt(2*gamma_pc/(gamma_pc - 1)*(1 - (Ps_Pc_pc)^((gamma_pc-1)/gamma_pc))) + 1*(Ps_Pc_pc - Pamb_Pc_pc);
C_E_cc = Gamma_cc*sqrt(2*gamma_cc/(gamma_cc - 1)*(1 - (Ps_Pc_cc)^((gamma_cc-1)/gamma_cc))) + eps_cc*(Ps_Pc_cc - Pamb_Pc_cc);

Isp_pc=c_est_pc*C_E_pc;
Isp_cc=c_est_cc*C_E_cc;
fun2 = subs(E == (mox_pc + mf_pc)*Isp_pc + (mox_cc + mf_cc)*Isp_cc);

%Criterios de Optimización


i = 1;
j = 1;
for OF_cc = OF_inicial:paso_OF:OF_final
    x_cc = OF_cc*M_H2/M_O2;
    OF_st = 8;
    if OF_cc > OF_st % Pobre -- Exceso de Oxidante
        a_cc = 1;
        b_cc = 0;
        c_cc = x_cc - 1/2;
        Mezcla = 'pobre';
    else % Rica
        a_cc = 2*x_cc;
        b_cc = 1-2*x_cc;
        c_cc = 0;
        Mezcla = 'rica';
    end
    for Pc_cc = Pc_inicial:paso_Pc:Pc_final
        syms mf_pc mf_cc
        funs = [subs(fun1), subs(fun2)];
        S = solve(funs);
        mf_cc = vpa(S.mf_cc);
        mf_pc = vpa(S.mf_pc);
        Isp(i,j) = double(vpa(subs(((mf_cc + mox_cc)*Isp_cc + mt*Isp_pc)/(mf + mox))));
        Km(i,j)=double(vpa(subs((K_ox*mox/mf+K_f)/(mox/mf+1))));
        rho_p(i,j)=double(vpa(subs((mox/mf+1)/((mox/mf)/densidad_ox+1/densidad_f))));
        %Delta_V(i,j)=double(vpa(subs(Isp(i,j)*log((1+Km(i,j)/(r+Km(i,j))))));
        % Gastos de oxidante
        c_estrella(i,j)=double(vpa(subs(c_est_cc)));
        C_empuje(i,j)=double(vpa(subs(C_E_cc)));
        
        MO(i,j) = double(vpa(subs(mox_cc + mox_pc)));
        MO_principal(i,j) = double(vpa(subs(mox_cc)));
        MO_sangrado(i,j) = double(vpa(subs(mox_pc)));
        
        % Gastos de oxidante
        MF(i,j) = double(vpa(subs(mf_cc + mf_pc)));
        MF_principal(i,j) = double(vpa(subs(mf_cc)));
        MF_sangrado(i,j) = double(vpa(subs(mf_pc)));
        
        j = j + 1;
    end
    i = i + 1;
    j = 1;
end
% end
vx = (meshgrid(Pc_inicial:paso_Pc:Pc_final)*1e-6);
vy = meshgrid(OF_inicial:paso_OF:OF_final)';


figure(1)
surf(vx,vy,Isp)
title('Isp')
xlabel('Pc [MPa]')
ylabel('O/F')
zlabel('Isp [m/s]')
colorbar

figure(2)
surf(vx,vy,Isp/(1+Km))
title('Isp/(1+Km)')
xlabel('Pc [MPa]')
ylabel('O/F')
zlabel('Isp/(1+Km) [m/s]')
colorbar

figure(3)
surf(vx,vy,Isp*rho_p)
title('Isp*rho_m')
xlabel('Pc [MPa]')
ylabel('O/F')
zlabel('Isp*Densidad [Kg/(m^2*s)]')
colorbar

%figure(4)
%surf(vx,vy,Isp*log((1+Km)/(r+Km)))
%title('Delta V')
%xlabel('Pc [MPa]')
%ylabel('O/F')
%zlabel('Delta V [m/s)]')
%colorbar

figure(5)
surf(vx,vy,c_estrella)
title('c*')
xlabel('Pc [MPa]')
ylabel('O/F')
zlabel('c* [m/s]')
colorbar

figure(6)
surf(vx,vy,C_empuje)
title('C_E')
xlabel('Pc [MPa]')
ylabel('O/F')
zlabel('C_E')
colorbar

figure(7)
surf(vx,vy,MO_principal)
title('Gasto Oxidante linea principal')
xlabel('Pc [MPa]')
ylabel('O/F')
zlabel('MO principal [kg/s]')
colorbar

figure(8)
surf(vx,vy,MO_sangrado)
title('Gasto Oxidante sangrado')
xlabel('Pc [MPa]')
ylabel('O/F')
zlabel('MO sangrado [kg/s]')
colorbar

figure(9)
surf(vx,vy,MF_sangrado)
title('Gasto Fuel sangrado')
xlabel('Pc [MPa]')
ylabel('O/F')
zlabel('MF sangrado [kg/s]')
colorbar

figure(10)
surf(vx,vy,MF_principal)
title('Gasto Fuel linea principal')
xlabel('Pc [MPa]')
ylabel('O/F')
zlabel('MF principal [kg/s]')
colorbar

% legend('Total','Principal','Sangrado')

% title('Isp/(1+km)')
% xlabel('O/F')
% ylabel('Pc [MPa]')
% zlabel('Isp [m/s]')
% colorbar