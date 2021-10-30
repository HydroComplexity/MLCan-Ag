function F = CAN_NIT_function( x )
%% Canopy Nitrogen distribution function

%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Rohit Nandan                                           %
%   Date         : April, 2020                                      %
%-------------------------------------------------------------------------%

% MLCan canopy model - conceptual hydrologic model

load ('NCAN_temp_variables.mat',...
            'PARAMS_temp', 'VARIABLES_temp', 'VERTSTRUC_temp', 'Cl', 'Cs',...
            'Cr', 'Cg', 'N_uptake_can_temp', 'CO2_Elev', 'si');

K_ec=x(1);
K_rub_leaf=x(2);
K_lhc_leaf=x(3);

LAI_ncan = VERTSTRUC_temp.LAIz;
Tl_sel=VARIABLES_temp.CANOPY.Tl(si);
Nrub_prev=sum(VARIABLES_temp.NCan.Nrub(si));
nl_can = PARAMS_temp.CanStruc.nl_can;

NstrF=PARAMS_temp.Nalloc.Nstr(si);
NlhcF=PARAMS_temp.Nalloc.Nlhc(si);
Rnchl=PARAMS_temp.Nalloc.Rnchl(si);


[Nleaf_tot_NU, Nr_NU, Ns_NU, Ng_NU, Nleaf_la_NU, Nrub_la_NU,...
    Vcmax_vert_opt, Chl_opt, fn] = VCHLcal(K_ec, K_rub_leaf, K_lhc_leaf,...
    N_uptake_can_temp, Nrub_prev,...
    LAI_ncan, Cl, Cr, Cs, Cg, Tl_sel, CO2_Elev, NstrF, NlhcF, Rnchl, nl_can);

VARIABLES_temp.NCan.Nall_start(si)=1;

VARIABLES_temp.CANOPY.PARabs_sun=VARIABLES_temp.NCan.Qabs_sun;
VARIABLES_temp.CANOPY.Tl_sun=VARIABLES_temp.NCan.Tl_sun;
VARIABLES_temp.CANOPY.Ci_sun=VARIABLES_temp.NCan.Ci_sun;

VARIABLES_temp.CANOPY.PARabs_shade=VARIABLES_temp.NCan.Qabs_shade;
VARIABLES_temp.CANOPY.Tl_shade=VARIABLES_temp.NCan.Tl_shade;
VARIABLES_temp.CANOPY.Ci_shade=VARIABLES_temp.NCan.Ci_shade;

Vcmax_vert_opt(isinf(Vcmax_vert_opt)|isnan(Vcmax_vert_opt)) = VARIABLES_temp.NCan.min_Vmax;

Chl_opt(isinf(Chl_opt)|isnan(Chl_opt)) = VARIABLES_temp.NCan.min_chl;

VARIABLES_temp.NCan.Vcmax_vert(:,si)=Vcmax_vert_opt;
VARIABLES_temp.NCan.Chl(si)=Chl_opt;

cntspecies=si;

sunlit=0;
[Ph_shade, An_shade ,Ph_limit, Jc_C4_shade, Jj_C4_shade, Js_C4_shade] = PHOTOSYNTHESIS_C4(VARIABLES_temp, PARAMS_temp, VERTSTRUC_temp, sunlit, cntspecies);

sunlit=1;
[Ph_sun, An_sun ,Ph_limit, Jc_C4_sun, Jj_C4_sun, Js_C4_sun] = PHOTOSYNTHESIS_C4(VARIABLES_temp, PARAMS_temp, VERTSTRUC_temp, sunlit, cntspecies);

LAI_sun =VARIABLES_temp.NCan.LAI_sun_can;
LAI_shade =VARIABLES_temp.NCan.LAI_shade_can;
fLAIz = VARIABLES_temp.NCan.fLAIz_can;

% An= sum(An_sun.*fLAIz,2).*LAI_sun + sum(An_shade.*fLAIz,2).*LAI_shade;
% F(1:nl_can) = -An(:);

Jj= sum(sum(Jj_C4_sun.*fLAIz,2).*LAI_sun + sum(Jj_C4_shade.*fLAIz,2).*LAI_shade);
Js= sum(sum(Js_C4_sun.*fLAIz,2).*LAI_sun + sum(Js_C4_shade.*fLAIz,2).*LAI_shade);
% Jc= sum(sum(Jc_C4_sun.*fLAIz,2).*LAI_sun + sum(Jc_C4_shade.*fLAIz,2).*LAI_shade);

F(1)=-Jj;
F(2)=-Js;
% F(3)=-Jc;

