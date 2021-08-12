%%%%%%%%%%%% MatlabToLatex %%%%%%%%%%%%%%%%%%%%%
clc, clear all, clear variables
%% PRECAMARA / TURBINA
% syms R_pc x_pc OF_pc Cp_o2_g Cp_h2_g Cp_h20_g M_h20 R_un M_o2 M_h2 M_prod_pc x_pc_sym Cp_prod_pc E_pc Cp_h20 C_pc Cp_h2 D_pc Cp_o2 E_pc C_pc D_pc T_c_prec T_ref Entalp_h2o E_pc C_pc D_pc E_pc Cp_prod_pc_sym Pi_t T_tin
%     latex(D_pc == x_pc_sym - 1/2)
%     latex(E_pc == 2*x_pc_sym)
%     latex(C_pc == 1-2*x_pc_sym)
% latex(M_prod_pc == (E_pc*M_h20 + C_pc*M_h2 + D_pc*M_o2)/(E_pc + C_pc + D_pc)) % g/mol
% latex(Cp_prod_pc == (E_pc*Cp_h20_g + C_pc*Cp_h2_g + D_pc*Cp_o2_g)/((E_pc + C_pc + D_pc))) % J/(g*K)
% latex(T_c_prec == T_ref - (Entalp_h2o*(E_pc/(C_pc+D_pc+E_pc))/Cp_prod_pc_sym)== T_tin)
% latex(OF_pc == x_pc * M_o2 / M_h2)
% latex(R_pc == R_un *1000/ M_prod_pc)
% latex(gamma_pc == Cp_prod_pc / (Cp_prod_pc - R_pc))

%% EC ACOPLAMIENTO
% syms Isp_0_aux Isp_0 Isp C_E_aux func_p_aux epsilon_aux ps_aux gamma_g_pc R_pc c_star_aux T_s_t Pi_t P_s_t Pi_pc Pi_iny_pc_ox P_e_t Gamma_gamma P_amb C_E P_s gamma_g epsilon c_star T_c x_cc Cp_medio_kg gamma_c Cp_medio_g Cp_medio R_c M_m C D E B A m_f_c m_o_c m_f_pc m_o_pc mt_m_v OF_pc OF_general P_df P_dox Pi_refr Pi_iny_cc_red Pi_iny_cc_ox Pc x mt_m m_bo_m Delta_P_bo rend_bo rho_ox  m_bf_m Delta_P_bf rend_bf rho_f rend_mec rend_t Cp_prod_pc T_tin gamma_pc
% 
% % latex(m_bf_m == 1 / (1 + OF_general))
% % latex(m_bo_m == OF_general/ (OF_general + 1))
% % latex(Delta_P_bo == (Pc / Pi_iny_cc_ox) - P_dox)
% % latex(Delta_P_bf == (Pc / (Pi_iny_cc_red * Pi_refr)) - P_df)
% 
% latex((((m_bo_m) * ( Delta_P_bo / (rend_bo * rho_ox) + m_bf_m * (Delta_P_bf / (rend_bf * rho_f))) == rend_mec * rend_t * mt_m * Cp_prod_pc * T_tin * ( 1 - (1 / Pi_t) ^ ( ( gamma_pc - 1) / gamma_pc )))))%, mt_m))
% 
% % latex(m_o_pc == mt_m_v * (OF_pc / (OF_pc + 1 )))
% % latex(m_f_pc == mt_m_v / (OF_pc + 1))
% % latex(m_o_c == m_bo_m - m_o_pc)
% % latex(m_f_c == m_bf_m - m_f_pc)
% % 
% % latex(A == m_f_c / 2)
% % latex(B == m_o_c / 32)
% % latex(E == min(A, 2*B))
% % 
% % latex(C == A - E) % if A > E
% % latex(C == 0) % else
% % 
% % latex(D == B - (1/2)* E) % if 2*B > E
% % latex(D == 0) % else
% % 
% % latex(M_m == ( E * M_h20 + (C * M_h2) + (D * M_o2)) / (C+D+E))
% % latex(R_c == R_un *1000/ M_m)
% % latex(Cp_medio == ( E * Cp_h20 + (C * Cp_h2) + (D * Cp_o2)) / (C+D+E))
% % latex( Cp_medio_g == (E*Cp_h20_g + C*Cp_h2_g + D*Cp_o2_g)/((E + C + D)))
% % latex(gamma_c == 1 / (1 - (R_c / Cp_medio_kg)))
% % latex(x_cc == (m_o_c/m_f_c) *M_o2 / M_h2)
% % latex(T_c == T_ref - (Entalp_h2o*(E/(C+D+E))/Cp_medio))
% % latex(c_star == sqrt(R_c * T_c))
% % func_p = (2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc)^((gamma_c-1)/gamma_c));
% % % latex(epsilon == (Gamma_gamma) / (((P_s/Pc) ^ (1 / gamma_c)) * sqrt((2*gamma_c/(gamma_c-1)) * (1 - (P_s/Pc)^((gamma_c-1)/gamma_c))))
% % latex(C_E == gamma_g* (sqrt(func_p)) + epsilon*((P_s/Pc) - (P_amb / Pc)))
% 
% %% TURBINA/TOBERA AUXILIAR
% % 
% % latex(P_e_t == (Pc * Pi_iny_pc_ox * Pi_pc) / ( Pi_iny_cc_ox ))
% % latex(P_s_t == P_e_t / Pi_t)
% % latex(T_s_t == T_tin * (1 - (rend_t*(1- Pi_t ^ ((gamma_pc-1)  / -gamma_pc))))) % SALE MAL
% 
% % latex(c_star_aux == sqrt(R_pc * T_s_t)/ gamma_g_pc)
% % 
% % latex(ps_aux == P_s_t / ( 1 + (gamma_pc-1)/2)^(gamma_pc/(gamma_pc-1)))
% % 
% % latex(epsilon_aux == gamma_g_pc / (((ps_aux/P_s_t) ^ (1 / gamma_pc)) * sqrt(func_p_aux)))
% % latex(C_E_aux == gamma_g_pc * sqrt(func_p_aux) + epsilon_aux*(((ps_aux/P_s_t)) - (P_amb / P_s_t)))
% % latex(Isp == Isp_0 * (1 - mt_m_v) + Isp_0_aux * mt_m_v)





%% OH NO HERE WE GO AGAIN
syms OF mf_C mo_C mf_PC mo_PC m_t mT_moPC moC_moPC Q Mp_cc OF_st OF_cc Cp_cc tau_bf tau_bo Cp_pc tau_t Mp_pc R_pc x_pc OF_pc Cp_o2_g Cp_h2_g Cp_h20_g M_h20 R_un M_o2 M_h2 M_prod_pc x_pc_sym Cp_prod_pc E_pc Cp_h20 C_pc Cp_h2 D_pc Cp_o2 E_pc C_pc D_pc T_c_prec T_ref Entalp_h2o E_pc C_pc D_pc E_pc Cp_prod_pc_sym Pi_t T_tin Isp_0_aux Isp_0 Isp C_E_aux func_p_aux epsilon_aux ps_aux gamma_g_pc R_pc c_star_aux T_s_t Pi_t P_s_t Pi_pc Pi_iny_pc_ox P_e_t Gamma_gamma P_amb C_E P_s gamma_g epsilon c_star T_c x_cc Cp_medio_kg gamma_c Cp_medio_g Cp_medio R_c M_m C D E B A m_f_c m_o_c m_f_pc m_o_pc mt_m_v OF_pc OF_general P_df P_dox Pi_refr Pi_iny_cc_red Pi_iny_cc_ox Pc x mt_m m_bo_m Delta_P_bo rend_bo rho_ox  m_bf_m Delta_P_bf rend_bf rho_f rend_mec rend_t Cp_prod_pc T_tin gamma_pc
%%%%%%%%% PRECAMARA %%%%%%%%%
% latex(Q ==  -(1 + 1/OF_st)*Entalp_h2o)
% 
% % %Pobre --> Exceso de Oxidante
% latex(Cp_pc == ( Cp_h20 + (OF_pc/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OF_pc/OF_st/2 - 1/2) ))
% latex(Mp_pc == ( M_h20 + (OF_pc/OF_st/2 - 1/2)*M_o2)/( 1 + (OF_pc/OF_st/2 - 1/2) ))
% latex(Cp_pc == Cp_pc/(Mp_pc*1e-3))
% latex((T_tin-T_ref) == (Q/Cp_pc)*(OF_st/(1+OF_pc)))
% % 
% % %Rica --> Exceso de Fuel
% latex(Cp_pc == ( Cp_h20 + (OF_st/OF_pc - 1)*Cp_h2 )/( 1 + (OF_st/OF_pc - 1) ))
% latex(Mp_pc == ( M_h20 + (OF_st/OF_pc - 1)*M_h2 )/( OF_st/OF_pc ))
% latex(Cp_pc == Cp_pc/(Mp_pc*1e-3))
% latex(T_tin == Q/Cp_pc*OF_pc/(1.+OF_pc) + T_ref)
% 
% latex(R_pc == (R_un / (Mp_pc*1e-3)))
% latex(gamma_pc == Cp_pc / (Cp_pc - R_pc))
% 
% %%%%%%%%% CAMARA COMBUSTION %%%%%%%%%
% latex(tau_t == rend_t*Cp_pc* T_tin*(1-(1/Pi_t)^((gamma_pc-1)/gamma_pc)))
% latex(tau_bo == Delta_P_bo / (rend_bo * rho_ox))
% latex(tau_bf == Delta_P_bf / (rend_bf * rho_f))
% %mezcla pobre
% latex(Cp_cc == ( Cp_h20 + (OF_cc/OF_st/2 - 1/2)*Cp_o2 )/( 1 + (OF_cc/OF_st/2 - 1/2)))
% latex(Mp_cc == ( M_h20 + (OF_cc/OF_st/2 - 1/2)*M_o2)/( 1 + (OF_cc/OF_st/2 - 1/2) ))
% latex(Cp_cc == Cp_cc/(Mp_cc*1e-3))
% latex(T_c == (Q/Cp_cc)*(OF_st/(1+OF_cc)) + T_ref)
% %mezcla rica
% latex(Cp_cc == ( Cp_h20 + (OF_st/OF_cc - 1)*Cp_h2 )/( 1 + (OF_st/OF_cc - 1) ))
% latex(Mp_cc == ( M_h20 + (OF_st/OF_cc - 1)*M_h2 )/( OF_st/OF_cc ))
% latex(Cp_cc == Cp_cc/(Mp_cc*1e-3))
% latex(T_c == Q/Cp_cc*OF_cc/(1+OF_cc) + T_ref)
% 
% latex(R_c == R_un*1000/ Mp_cc)
% latex(gamma_c == 1 / (1 - (R_c / Cp_cc)))
% 
% %%%%%%%%% TURBINA/TOBERA AUXILIAR %%%%%%%%%
% %%%%%%%%% ACOPLAMIENTO MECANICO %%%%%%%%%
% latex(moC_moPC == ( (1 + 1/OF_pc)*rend_mec* tau_t - (tau_bo + tau_bf/OF_pc) )/...
%                (tau_bo + tau_bf/OF_cc))
% latex(mT_moPC == moC_moPC*(1 + 1/OF_cc) + (1 + 1/OF_pc))
% latex(Isp == Isp_0 * (1 -(1 + 1/OF_pc)/mT_moPC) + Isp_0_aux * (1 + 1/OF_pc)/mT_moPC)

latex( m_t == E / Isp )
latex( mo_PC == m_t / mT_moPC )
latex( mf_PC == mo_PC / OF_pc )
latex( mo_C == mo_PC * moC_moPC )
latex( mf_C == mo_C / OF_cc)
latex( OF == (mo_C + mo_PC)./(mf_C + mf_PC))

