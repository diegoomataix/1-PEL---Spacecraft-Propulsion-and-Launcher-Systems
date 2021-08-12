clear all, clear variables, close all

%% DATOS
M_cp = 100;          % kg
M_so = 20;           % kg
Mppo = 1;            % kg

delta = 0.03; 
alpha_pp = 20 / 10^3;       % kg / kW --> % kg / W
alpha_m = 4 / 10^3;         % kg / kW --> % kg / W

el = 1.602176565e-17; %100 * 1.60218E-19;   % 100 eV (pérdidas por ion) --> J
eff_0 = 0.75;
N_A = 6.022E23;
qi = 1.6E-19;        % C
perm_vac = 8.81E-12; % C^2 N^-1 m^-2

M_Cs = 132.9;        % g/mol
M_Xe = 131.3;
M_Hg = 200.6;

m_Cs = (M_Cs*10^(-3)) / N_A;    %2.207E-25;    % kg
m_Xe = (M_Xe*10^(-3)) / N_A;    %2.180E-25;    % kg
m_Hg = (M_Hg*10^(-3)) / N_A;    %3.331E-25;    % kg

I_tot = 9E5;         % Ns

coste_kg_A = 50000;
coste_kg_B = 9000;      % $/kg
Delta_V = 4760;      % m/s

M_0a = 100 + 20 + 1;

L = 0.7 * 10^-3;     % m
A = 30 * 10^-4;

%% B)OPT
coste_A = M_0a * coste_kg_A;

Isp = linspace(1000*9.81,3000*9.81,20);

for i = 1:length(Isp)
     
M_p(i) = I_tot / Isp(i);
M_0(i) = M_p(i) / ( 1  - exp(-Delta_V/Isp(i)));
P(i) = (M_0(i) - M_0a - (1 + delta) * M_p(i)) / (alpha_pp + alpha_m);
coste_B(i) = coste_kg_B * M_0(i);

% Cesio
eff_Cs(i) = eff_0 * (( Isp(i)^2 ) / (( Isp(i)^2) + 2 *( el / m_Cs)));
E_P_Cs(i) = 2 * eff_Cs(i) / Isp(i);
E_Cs(i) = E_P_Cs(i) * P(i);
tb_Cs(i) = I_tot / E_Cs(i);
tb_dias_Cs(i) = tb_Cs(i) / (60 * 60 * 24);

% Xenon
eff_Xe(i) = eff_0 * ( Isp(i)^2 ) / ( Isp(i)^2 + 2 *( el / m_Xe));
E_P_Xe(i) = 2 * eff_Xe(i) / Isp(i);
E_Xe(i) = E_P_Xe(i) * P(i);
tb_Xe(i) = I_tot / E_Xe(i);
tb_dias_Xe(i) = tb_Xe(i) / (60 * 60 * 24);

% Mercurio
eff_Hg(i) = eff_0 * ( Isp(i)^2 ) / ( Isp(i)^2 + 2 *( el / m_Hg));
E_P_Hg(i) = 2 * eff_Hg(i) / Isp(i);
E_Hg(i) = E_P_Hg(i) * P(i);
tb_Hg(i) = I_tot / E_Hg(i);
tb_dias_Hg(i) = tb_Hg(i) / (60 * 60 * 24);

end

%% PLOT

[i, j] = find(ismember(coste_B, min(coste_B(:))))   
% 
figure()
plot(Isp,P,'LineWidth',2)
hold on
xlabel('Isp [m/s]')
ylabel('P [W]')
grid on
hold off

figure()
plot(Isp,M_0,'LineWidth',2)
hold on
xlabel('Isp [m/s]')
ylabel('M_0 [kg]')
grid on
hold off

figure()
plot(Isp,(coste_B/10^6),'LineWidth',2)
hold on
xlabel('Isp [m/s]')
ylabel('Coste [M$]')
grid on
hold off

figure()
plot(Isp,(E_P_Cs*10^6),'LineWidth',2,'Color','red')
hold on
plot(Isp,(E_P_Xe*10^6),'LineWidth',2, 'Color','green')
plot(Isp,(E_P_Hg*10^6),'LineWidth',2, 'Color','blue')
grid on
legend('Cs','Xe','Hg')
xlabel('Isp [m/s]')
ylabel('E/P [mN/kW]')
hold off

figure()
plot(Isp,tb_dias_Cs,'LineWidth',2,'Color','red')
hold on
plot(Isp,tb_dias_Xe,'LineWidth',2,'Color','green')
plot(Isp,tb_dias_Hg,'LineWidth',2,'Color','blue')
legend('Cs','Xe','Hg')
xlabel('Isp [m/s]')
ylabel('tb [dias]')
grid on
hold off

figure()
plot(Isp,eff_Cs,'LineWidth',2,'Color','red')
hold on
plot(Isp,eff_Xe,'LineWidth',2,'Color','green')
plot(Isp,eff_Hg,'LineWidth',2,'Color','blue')
legend('Cs','Xe','Hg')
xlabel('Isp [m/s]')
ylabel('\eta_m')
grid on
hold off

%% 
%nos falta sacar:
    % Mpp
    % I
    % Ms
    % gasto másisco
    
%% C) 3 ELECTRODOS 
Isp_opt = max(Isp(:));

Va = m_Hg * Isp_opt^2 / (2 * qi);

E_A = E_Hg(20) / A; 
syms Delta_Volt_sym
Delta_Volt = vpasolve(E_A == (8/9) * perm_vac * ((Delta_Volt_sym / L)^2) * sqrt(Va/Delta_Volt_sym));


%% BACKUP
% % m = [m_Cs, m_Xe, m_Hg]
% % 
% % for h = 1:length(m)
%    
% Isp = linspace(1000*9.81,3000*9.81,20);
% %tb = linspace((7776000/3),31104000, 20);      % 3 meses - 12 meses
% 
% for i = 1:length(Isp)
%      
% M_p(i) = I_tot / Isp(i);
% M_0(i) = M_p(i) / ( 1  - exp(-Delta_V/Isp(i)));
% P(i) = (M_0(i) - M_0a - (1 + delta) * M_p(i)) / (alpha_pp + alpha_m);
% 
% eff(i) = eff_0 * ( Isp(i)^2 ) / ( Isp(i)^2 + 2 *( el / m_Cs));
% E_P(i) = 2 * eff(i) / Isp(i);
% E(i) = E_P(i) * P(i);
% tb(i) = I_tot / E(i);
% tb_dias(i) = tb(i) / (60 * 60 * 24);
% 
% % for j = 1: length(tb)
% %     
% % eff(i) = eff_0 * ( Isp(i)^2 ) / ( Isp(i)^2 + 2 *( el / m_Cs));
% % 
% % E_P(i) = 2 * eff(i) / Isp(i);
% % 
% % P_chi_E(i) = Isp(i) / 2;
% % 
% % I_E(i) = (qi / m_Cs) * (1 / Isp(i));
% % 
% % m_dot_E(i) = 1 / Isp(i);
% %     
% % E(j) = I_tot / tb(j);
% % 
% % P(i, j) = E(j) / E_P(i); 
% % 
% % M_0b(i,j) = (M_0a + P(i,j) * ( alpha_pp + alpha_m)) / ( 1 - ( 1+ delta) * (1 - exp(- Delta_V/ Isp(i))));
% % 
% % coste_B(i,j) = coste_kg_B * M_0b(i,j);
% % 
% %     end
% end

