clc;
clear all;
close all;


%% DATOS

% Rendimientos
eta_m = 0.94;
eta_t = 0.68;
etab_H = 0.62; etab_O = 0.64;

% Presiones
Pd_H = 3e5;
Pd_O = 6e5;

piiny_ccO = 0.78;
piiny_ccH = 0.82;
piref = 0.88;
piiny_pcO = 0.74;
pipc = 0.92;

% Props quimicas
ro_H = 71;
ro_O = 1140;

MO = 15.999; MO2 = MO*2;
MH = 1.008; MH2 = MH*2;
MH2O = MH2 + MO;

cpO2 = 39.07/(MO2*1e-3); 
cpH2 = 36.11/(MH2*1e-3); 
cpH2O = 54.32/(MH2O*1e-3);

R = 8314.4598;
gammaO2 = cpO2/(cpO2 - R/MO2);
gammaH2 = cpH2/(cpH2 - R/MH2);
gammaH2O = cpH2O/(cpH2O - R/MH2O);

% Combustion
hfH2O = -241818e3/(MH*2 + MO);
Tref = 298.15;

% Turbina
pit_max = 12;
Te_t = 900;

% Masas
ko = 0.02;
kf = 0.2;
r = 0.2;

% Actuaciones
E = 50e3*9.81;
Ps = 40e3;
Pa = 101325;
%Pa = Ps;


%% TEMPERATURA DE COMBUSTIÓN EN FUNCION DE OF

% Relaciones de mezcla
OFst = 8;
cpH2 = 36.11; cpO2 = 39.07; cpH2O = 54.32;  % mol

% Combustion en funcion de OF
[OFc, cp_prod, M_prod, T_prod, gamma_prod, Gamma_prod, cest_prod ] = ...
         Combustion(OFst, cpH2, cpO2, cpH2O, MH2, MO2, MH2O, R, hfH2O, Tref);



     
%% PRE-COMBUSTOR 
[minimo, pos_t] = min( abs(T_prod-Te_t) );
OFpc = OFc(pos_t);
cp_t = cp_prod(pos_t);
gamma_t = gamma_prod(pos_t);
Gamma_t = Gamma_prod(pos_t);

% figure();
%     hold on
%     plot(OFc, T_prod,'LineWidth',2)
%     plot([0,16], [900,900],'k','LineWidth',2)
%     title(['Relacion OF_{PC} = ' num2str(OFpc) ' para T_{et} = 900 K'])
%     xlabel('OF precamara')
%     ylabel('T_{et}')
%     legend('T_C', '900K')
%     grid on
%     box on


%% TURBOMAQUINARIA
% Bombas
BO = Bomba(etab_O, ro_O, Pd_O, piiny_ccO);
BH = Bomba(etab_H, ro_H, Pd_H, piiny_ccH*piref);

% Turbina
Turbina = Turbina(eta_t, pit_max, Te_t, gamma_t, cp_t);

%% TOBERA CONVERGENTE

epsilon_t = 1; 
cest_t = cest(R, M_prod(pos_t), Turbina.Ts, Gamma_prod(pos_t));
pstps = ( (gamma_t+1)/2 )^( gamma_t/(gamma_t-1) ); % pstps = ( (gamma_t+1)/2 )^( gamma_t/(gamma_t-1) );                %relación entre presió a la entrada tobera( salida turbina ) y presion salida tobera (condición de bloqueo)
pspst = 1/pstps;


%% BARRIDO PRESIONES

Pc =  [100e5:1e5:700e5];%32.892476699999996e6;

%OFc = input('OF de la Camara principal = ');
% Combustion en funcion de OF
[cp_prod, M_prod, Tc_prod, gamma_prod, Gamma_prod, cest_prod ] = ...
         Combustion_OF_Fijo(OFst, OFc, cpH2, cpO2, cpH2O, MH2, MO2, MH2O, R, hfH2O, Tref);


%Pc = input('Presion de la Camara principal = ');
for i = 1:length(Pc)
     
    % Relacion de areas funcion de relacion de presiones
    epsilon = Epsilon( gamma_prod, Gamma_prod, Ps/Pc(i) );
    
    % Isp tobera principal ADAPTADO
    cEo = cE(gamma_prod, Gamma_prod, epsilon, Ps/Pc(i),Pc(i), Pa);
    Ispo = cest_prod.*cEo;   
    
    % TURBOMAQUINARIA
    BO.Punto_Funcionamiento(Pc(i));
    BH.Punto_Funcionamiento(Pc(i));
    
    % Turbina
    Ps_pc = BO.Ps*piiny_pcO*pipc;                     %presión a la salida del postcombustor (entrada turbin
    Turbina.Punto_Funcionamiento(Ps_pc)                                          %presion a la salida de la turbina
    cE_t = cE(gamma_t, Gamma_t, epsilon_t, pspst, Turbina.Ps, Pa);
    Ispt = cest_t*cE_t;                                               %Impulso turbina
    
    % Acoplamiento turbomaquinaria
    moC_moPC = ( (1 + 1/OFpc)*eta_m*Turbina.tau - (BO.tau + BH.tau./OFpc) )./...
               (BO.tau + BH.tau./OFc);
    mT_moPC = moC_moPC.*(1 + 1./OFc) + (1 + 1/OFpc);

    % Isp total
    Isp_T(i,:) = ( 1 - (1 + 1/OFpc)./mT_moPC ).*Ispo + (1 + 1/OFpc)./mT_moPC*Ispt;

 
end

%% GASTO MASICO TOTAL

mT = E./Isp_T;

moPC = mT./mT_moPC; 
mfPC = moPC/OFpc;

moC = moPC.*moC_moPC; 
mfC = moC/OFc;

OFT = (moC + moPC)./(mfC + mfPC);

% RESULTADOS


Resultados(1).OF_Global = OFT;
Resultados(1).OF_CC = OFc;
Resultados(1).OF_PC = OFpc;
Resultados(1).P_amb = Pa;
Resultados(1).Pc = Pc;
Resultados(1).Tc = Tc_prod;
Resultados(1).cp_prod = cp_prod;
Resultados(1).M_prod = M_prod;
Resultados(1).Tc_prod = Tc_prod;
Resultados(1).gamma_prod = gamma_prod;
Resultados(1).Gamma_prod = Gamma_prod;
Resultados(1).cest_prod = cest_prod;

Tabla_Resultados = struct2table(Resultados)



% APARTADO 1
{
Isp_T_max == max(Isp_T, [], 'All');
[fila, columna] == find(Isp_T == Isp_T_max);

Pc_opt = Pc(fila);
OF_opt = OFc(columna);

figure()
    s = surf(OFc,Pc,Isp_T);
    s.EdgeColor = 'none';
    xlabel('OF')
    ylabel('Pc')
    zlabel('Isp')
    
figure();
    hold on
    for i = 1:length(Pc)
        plot(OF, Isp_T(i,:))
        [v, p] = max(Isp_T(i,:));
        title(['max = ' num2str(v), 'OF = ' num2str(OF(p))])
    end
    legend(num2str(Pc'/1e5))
    xlabel('OF')
    ylabel('Isp T')
    title('Isp OF para varias Pcs')
    
    disp(['Isp max = ' num2str(Isp_T_max)])
    disp(['OF para Isp max = ' num2str(OF_opt)])
    disp(['Pc para Isp max = ' num2str(Pc_opt)])
    
   
    
% APARTADO 2


k = (ko*OFT + kf)./(OFT + 1);
criterio_2 = Isp_T./(1 + k);    % es un + o -?

criterio_2_max = max(criterio_2, [], 'All');
[fila, columna] = find(criterio_2 == criterio_2_max);

disp(['OF para Isp max = ' num2str(OFc(columna))])
disp(['Pc tal que Isp max = ' num2str(Pc(fila))])

figure()
    s = surf(OFT,Pc,criterio_2);
    s.EdgeColor = 'none';
    xlabel('OF')
    ylabel('Pc')
    zlabel('criterio_2')
    title('CRITERIO 2')

% APARTADO 3

ro = (OFT.*ro_O + ro_H)./(1 + OFT); % (OFT + 1)./( OFT./ro_H + 1/ro_O );
Impulso_volumetrico = Isp_T.*ro;

Impulso_volumetrico_max = max(Impulso_volumetrico, [], 'All');
[fila, columna] = find(Impulso_volumetrico == Impulso_volumetrico_max);

disp(['Isp volumetrico max = ' num2str(Impulso_volumetrico_max)])
disp(['OF para Isp max = ' num2str(OFc(columna))])
disp(['Pc tal que Isp max = ' num2str(Pc(fila))])

figure()
    s = surf(OFT,Pc,Impulso_volumetrico);
    s.EdgeColor = 'none';
    xlabel('OF')
    ylabel('Pc')
    zlabel('criterio_3')
    title('Isp DEPOSITIOS')
 

% APARTADO 4

delta_v = Isp_T.*log( (1 + k)./(r + k) );

delta_v_max = max(delta_v, [], 'All');
[fila, columna] = find(delta_v == delta_v_max);

disp(['OF para Isp max = ' num2str(OFc(columna))])
disp(['Pc tal que Isp max = ' num2str(Pc(fila))])

figure()
    s = surf(OFT,Pc,delta_v);
    s.EdgeColor = 'none';
    xlabel('OF')
    ylabel('Pc')
    zlabel('criterio_4')
    title('DELTA V')
}