clear all, clear variables

%b
 syms L A E_A Isp_opt qi perm_vac Va Delta_Volt m_Hg tb P Delta_V E_P E M_p I_tot Isp M_0 M_0a delta el m alpha_m alpha_pp coste_kg_B coste_B eff_0 eff
% 
% latex(M_p == I_tot / Isp)
% 
% latex(M_0 == M_p / ( 1  - exp(-Delta_V/Isp)))
% 
% latex(P == (M_0 - M_0a - (1 + delta) * M_p) / (alpha_pp + alpha_m))
% 
% latex(coste_B == coste_kg_B * M_0)
% 
% latex(eff == eff_0 * (( Isp^2 ) / (( Isp^2) + 2 *( el / m))))
% 
% latex(E_P == 2 * eff / Isp)
% 
% latex(E == E_P * P)
% 
% latex(tb == I_tot / E)


%c

latex(Va == m * Isp_opt^2 / (2 * qi))

latex(E_A == E / A)

latex(E_A == (8/9) * perm_vac * ((Delta_Volt / L)^2) * sqrt(Va/Delta_Volt))