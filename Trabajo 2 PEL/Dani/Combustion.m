function [OF, cp_prod, M_prod, Tc_prod, gamma_prod, Gamma_prod, cest ] = ...
         Combustion(OFst, cpH2, cpO2, cpH2O, MH2, MO2, MH2O, R, hfH2O, Tref)


    % Relaciones de mezcla
    OFmin = 0.5;
    OFmax = 16;
    OFp = linspace(OFmin,OFst,1e3);
    OFr = linspace(OFst,OFmax,50);

    % Calculo de temperatura de combustion
    Q = -(1 + 1/OFst)*hfH2O;

    % cp promedio
    OFr = 0.7477;
    cpp = ( cpH2O + (OFst./OFp - 1)*cpH2 )./( 1 + (OFst./OFp - 1) );
    cpr = ( cpH2O + (OFr/OFst/2 - 1/2)*cpO2 )./( 1 + (OFr/OFst/2 - 1/2) );

    % Masa molar promedio
    Mpp = ( MH2O + (OFst./OFp - 1)*MH2 )./( OFst./OFp );
    Mpr = ( MH2O + (OFr/OFst/2 - 1/2)*MO2)./( 1 + (OFr/OFst/2 - 1/2) );

    % cp promedio en kg
    cpp = cpp./(Mpp*1e-3);  
    cpr = cpr./(Mpr*1e-3);

    % Temperatura de combustion
    Tcp = Q./cpp.*OFp./(1.+OFp) + Tref;
    Tcr = Q./cpr.*OFst./(1.+OFr) + Tref;

    % gamma promedio en kg
    gamma_p = cpp./(cpp - R./Mpp);
    gamma_r = cpr./(cpr - R./Mpr);

    % Propiedades de los productos
    OF = [OFp, OFr];
    cp_prod = [cpp cpr];
    M_prod = [Mpp, Mpr];
    Tc_prod = [Tcp, Tcr];
    gamma_prod = [gamma_p gamma_r];
    Gamma_prod = Gamma(gamma_prod);

    % c* en funcion OF
    cest = sqrt(R./M_prod.*Tc_prod)./Gamma_prod;

    % Entalpias
    %{
    Tcr = -hfH2O./( cpH2O + MO2/MH2O.*(OFr/OFst/2 - 1/2)*cpO2 ) + Tref;
    Tcp = -hfH2O./( cpH2O + MH2/MH2O.*(OFst./OFp - 1)*cpH2 ) + Tref;

    figure();
        hold on
        plot(OFp, Tcp,'LineWidth',2)
        plot(OFr, Tcr,'LineWidth',2)
        plot([0,16], [900,900],'k','LineWidth',2)
        title('T ENTRADA TURBINA vs OF')
        xlabel('OF precamara')
        ylabel('T_{et}')
        legend('Pobre', 'Rica', '900K')
        grid on
        box on
    %}    


end

