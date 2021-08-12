function Gamma = Gamma(gamma)

    Gamma = sqrt(gamma).*(2./(gamma+1)).^((gamma+1)./2./(gamma-1));
    
end

