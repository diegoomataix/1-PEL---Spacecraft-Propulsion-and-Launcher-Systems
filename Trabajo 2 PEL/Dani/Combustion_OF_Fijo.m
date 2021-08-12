function [cp_prod, M_prod, Tc_prod, gamma_prod, Gamma_prod, cest ] = ...
         Combustion_OF_Fijo(OFst, OF, cpH2, cpO2, cpH2O, MH2, MO2, MH2O, R, hfH2O, Tref)


    % Relaciones de mezcla
    OFmin = 0.5;
    OFmax = 16;
    OFp = OF;


    % Calculo de temperatura de combustion
    Q = -(1 + 1/OFst)*hfH2O;

    % cp promedio
    cpp = ( cpH2O + (OFst./OFp - 1)*cpH2 )./( 1 + (OFst./OFp - 1) );

    % Masa molar promedio
    Mpp = ( MH2O + (OFst./OFp - 1)*MH2 )./( OFst./OFp );

    % cp promedio en kg
    cpp = cpp./(Mpp*1e-3);  

    % Temperatura de combustion
%     Tcp = Q./cpp.*OFp./(1.+OFp) + Tref;
    Tcp = Q./cpp.*OFp./(1.+OFp) + Tref;

    % gamma promedio en kg
    gamma_p = cpp./(cpp - R./Mpp);

    % Propiedades de los productos
    OF = OFp;
    cp_prod = cpp;
    M_prod = Mpp;
    Tc_prod = [Tcp];
    gamma_prod = [gamma_p];
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

