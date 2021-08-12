function cE = cE(gamma, Gamma, epsilon, pspct, pc, pa)

    cE = Gamma.*sqrt(2.*gamma./(gamma-1).*(1 - (pspct).^((gamma-1)./gamma))) + ...
          epsilon.*(pspct - pa./pc);

end

