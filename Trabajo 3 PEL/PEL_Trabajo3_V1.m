%% DATOS

masa_cp = 1500;             %[kg]
masa_prop = 118;            %[kg]

emp_1 = 8;                  %[N], motor principal
emp_2 = 0.7;                %[N], motor secundario

x = 0.4;

rho_N2H4 = 1032;            %[kg/m^3]
P_c = 10*10^5;              %[Pa]
Pi_tobera = 80;             %Pc/Ps

DeltaP_iny = 0.1;

Cp_H2 = 29.9;               %[J/mol*K]
Cp_N2 = 32.1;               %[J/mol*K]
Cp_NH3 = 54.0;              %[J/mol*K]
Cp_N2H4 = 98.8;             %[J/mol*K]

Hf_N2H4 = 50.6*10^3;        %[J/mol]
Hf_NH3 = -45.9*10^3;        %[J/mol]

M_N = 14*10^(-3);           %[kg/mol]
M_H = 1*10^(-3);            %[kg/mol]

T_ref = 298.15;             %[K]
T_iny = T_ref;
R_u = 8.31446;              %[J/mol*K]


%% Apartado a

%ecuacion de combustion y decomposicion del N2H4
A = 1;                                                              %N2H4
B = (4*(1-x))/3;                                                    %NH3
C = (2*x+1)/3;                                                      %N2
D = 2*x;                                                            %H2

M_prod = (B*(M_N + 3*M_H) + C*2*M_N + D*2*M_H)/(B+C+D);                     %Masa molecular de los productos [kg/mol]
Cp_p = (B*Cp_NH3 + C*Cp_N2 + D*Cp_H2)/(B+C+D);                              %Cp medio [J/mol*K]
Cp_p_m = Cp_p / M_prod;                                                     %Cp medio [J/kg*K]
R_p = R_u / M_prod;                                                         %R productos [J/kg*K]
Gamma = 1/(1-(R_p/Cp_p_m));                                                 
G_Gamma = sqrt(Gamma)*(2/(Gamma+1))^((Gamma+1)/(2*(Gamma-1)));

DeltaT = (Hf_N2H4 - B*Hf_NH3)/(B*Cp_NH3 + C*Cp_N2 + D*Cp_H2);
T_c = DeltaT + T_ref;

c_star = sqrt(R_p*T_c)/G_Gamma;

%Tobera
epsilon = G_Gamma / (((1/Pi_tobera)^(1/Gamma))*sqrt((2*Gamma/(Gamma-1))*(1-((1/Pi_tobera)^((Gamma-1)/Gamma)))));
%satelite de comunicacion  : actua en el vacio
Ce = G_Gamma*sqrt((2*Gamma/(Gamma-1))*(1-((1/Pi_tobera)^((Gamma-1)/Gamma)))) + epsilon*(1/Pi_tobera);

%Gastos masico y Areas de garganta
Ag_1 = emp_1/(P_c*Ce);                                                      %Area de garganta motor principal [m^2]
Ag_2 = emp_2/(P_c*Ce);                                                      %Area de garganta motor segundario [m^2]

m_1 = (P_c*Ag_1)/c_star;                                                    %gasto masico motor principal [kg/s]
m_2 = (P_c*Ag_2)/c_star;                                                    %gasto masico motor segundario [kg/s]


%% Apartado b
P_t = P_c + DeltaP_iny*P_c;                                                 %presion en el tanque [Pa]
V_t = masa_prop / rho_N2H4;                                                 %volumen del tanque [m^3]


%% Apartado c
m_out = m_1 + m_2;                                                          %gasto masico que sale del tanque
alpha = 0.15;                                                               

%Propriedades de gas de presurizacion y del deposito (diferentes configuraciones)
%(Aluminio, Acero, Material Compuesto)
rho_dep = [2810 8190 1600];                                                 %[kg/m^3]
sigma_dep = [0.45*10^9 0.540*10^9 3.8*10^9];                                %[Pa]

%(H2, He, Aire, CO2, Ar)
M_gas = [2.014*10^(-3) 4.002*10^(-3) 28.97*10^(-3) 44.01*10^(-3) 39.94*10^(-3)];     %[kg/mol]
Gamma_gas = [1.41 1.66 1.4 1.3 1.67];
R_gas = [0 0 0 0 0];                                                        %initializacion R_gas

T_d0 = 294;                                                                 %temperatura del deposito al instante initial
%!! si cambia la temperatura cambian todas las propriedades
N = 15;                                                                     %numero de valores de P_d0
P_d0_min = (1+alpha)*P_t;                                                   %vallor minimal en el linspace de P_d0
P_d0 = linspace(13*10^5,27*10^5,N);                                        %Presion del deposito al instante initial, P_d0 min = (1+alpha+Gamma_gas*V_tt/vd)*P_t

gases = 5;                                                                  %gas de presurizacion considerado (H2, He, Aire, CO2, Ar)
masa_gas = zeros(gases, N);
Vd = zeros(gases, N);

masa_alu = zeros(gases, N);
masa_t_alu = zeros(gases, N);
masa_acero = zeros(gases, N);
masa_t_acero = zeros(gases, N);
masa_mc = zeros(gases, N);
masa_t_mc = zeros(gases, N);

for j = 1:N
    for i = 1:gases
        R_gas(i) = R_u/M_gas(i);
        
        Vd(i,j) = Gamma_gas(i)*P_t*V_t/(P_d0(j)*(1-(1+alpha)*P_t/P_d0(j)));
        
        masa_gas(i,j) = P_d0(j)*Vd(i,j)/(R_gas(i)*T_d0);                    %masa del gaz en funcion del gaz y de P_d0
        
        
        masa_alu(i,j) = 3*P_d0(j)*Vd(i,j)*rho_dep(1)/(2*sigma_dep(1));
        masa_t_alu(i,j) = masa_gas(i,j) + masa_alu(i,j);
        
        masa_acero(i,j) = 3*P_d0(j)*Vd(i,j)*rho_dep(2)/(2*sigma_dep(2));
        masa_t_acero(i,j) = masa_gas(i,j) + masa_acero(i,j);
        
        masa_mc(i,j) = 3*P_d0(j)*Vd(i,j)*rho_dep(3)/(2*sigma_dep(3));
        masa_t_mc(i,j) = masa_gas(i,j) + masa_mc(i,j);
        
    end
end











