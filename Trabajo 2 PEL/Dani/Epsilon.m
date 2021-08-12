function Epsilon = Epsilon(gamma_garg, Gamma_garg, pspcs)

    Epsilon = Gamma_garg./( ( pspcs.^(1./gamma_garg) ).*...
        sqrt(2.*gamma_garg./(gamma_garg-1).*...
        (1 - (pspcs).^((gamma_garg-1)./gamma_garg))) ) ;

end

