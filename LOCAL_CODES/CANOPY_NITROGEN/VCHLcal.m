function [Nleaf_tot, Nr, Ns, Ng, Nleaf_la, Nrub_la, Vcmax_vert_opt, Chl_opt, fn] =VCHLcal(K_ec, K_rub_leaf, K_lhc_leaf, N_uptake_can, Nrub_prev,...
    LAI_ncan, Cl, Cr, Cs, Cg, Tl_sel, CO2_Elev, NstrF, NlhcF, Rnchl, nl_can)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

min_Vmax=8;
max_Vmax=65;

N_all=N_uptake_can;

gmtommol=1/0.0140067;
mmoltogm=0.0140067;

LAD=...
    [0.0369127665665577; 0.0395397884964049; 0.0420898857612072;...
    0.0445252956223973; 0.0468081563454242; 0.0489014705498921;...
    0.0507700933998415; 0.0523817091634002; 0.0537077584441472;...
    0.0547242790146216; 0.0554126257063032; 0.0557600391530863;...
    0.0557600391530863; 0.0554126257063032; 0.0547242790146216;...
    0.0537077584441472; 0.0523817091634002; 0.0507700933998415;...
    0.0489014705498921; 0.0468081563454242];

%% LAI cum, Biomass and fraction of C allocation
LAI_ncan_cum = cumsum ( flipdim(LAI_ncan(:),1) );
LAI_ncan_cum_top = flipdim(LAI_ncan_cum(:), 1);
LAI_ncan_sum=sum(LAI_ncan);

Fleaf=Cl/(Cl+Cs+Cr+Cg); Fleaf(isinf(Fleaf)|isnan(Fleaf)) = 0;
Froot=Cr/(Cl+Cs+Cr+Cg); Froot(isinf(Froot)|isnan(Froot)) = 0;
Fstem=Cs/(Cl+Cs+Cr+Cg); Fstem(isinf(Fstem)|isnan(Fstem)) = 0;
Fgrain=Cg/(Cl+Cs+Cr+Cg); Fgrain(isinf(Fgrain)|isnan(Fgrain)) = 0;

%% Nitrogen Allocation

if(Cg>0)
    fn=0;
    f_g=0.1 * exp(2*Fgrain);
    Nr = Froot * N_all * (1-f_g);
    Ns = Fstem * N_all * (1-f_g);    
    Nleaf_tot = Fleaf * N_all * (1-f_g);
    Ng = (Fgrain * N_all)+((Fleaf + Froot+ Fstem)*N_all)*(f_g);
    
else
    fn = 0.2 * exp(-2*Nrub_prev);
    Nr = Froot * (1-fn) * N_all;
    Ns = Fstem * (1-fn) * N_all;
    Nleaf_tot = (Fleaf * N_all) + (fn * (Froot * N_all + Fstem * N_all));
    Ng = Fgrain * N_all;
    
end

% Nleaf_tot = (Fleaf * N_all) + (fn * (Froot * N_all + Fstem * N_all));
% Ng = Fgrain * N_all;

Nleaf_tot(isnan(Nleaf_tot)|isnan(Nleaf_tot))=0;
Nleaf_tot(Nleaf_tot<0)=0;

Nleaf_tot_mmol = Nleaf_tot * gmtommol;

% Nb_mmol=12.5;       %Anten et al 1995 Structural Nitrogen mmol N/ga m2

Nb_mmol=NstrF * Nleaf_tot_mmol* (1-fn)/LAI_ncan_sum;

% Hirose 2005
Nleaf_mmol=(((K_ec*(Nleaf_tot_mmol-(Nb_mmol*LAI_ncan_sum))*exp(-K_ec*LAI_ncan_cum_top))/(1-exp(-K_ec*LAI_ncan_sum)))+Nb_mmol)/nl_can;

% uniform solution
% Nleaf_mmol=(Nleaf_tot_mmol/nl_can)*ones(nl_can,1);

Nleaf_la=Nleaf_mmol*mmoltogm;

%% Vmax
Nrub_la_mmol = (K_rub_leaf * Nleaf_mmol) + (fn * (Froot + Fstem)* N_all + fn * Nb_mmol*mmoltogm ).*(LAD./LAI_ncan_sum);

Nrub_ga_mmol=Nrub_la_mmol.*LAI_ncan;

% Vcmax_vert_opt=(b*Nrub_la_mmol.*Temp_fact)./(1+(0.04*Temp_fact.^3));
if CO2_Elev==1
    Vcmax_vert_opt=(6.25*Nrub_la_mmol*27.5*0.014)*20;
else
    Vcmax_vert_opt=(6.25*Nrub_la_mmol*20.7*0.014)*20;
end

Vcmax_vert_opt(isinf(Vcmax_vert_opt)|isnan(Vcmax_vert_opt)) = min_Vmax;
Vcmax_vert_opt(Vcmax_vert_opt<min_Vmax)=min_Vmax;
Vcmax_vert_opt(Vcmax_vert_opt>max_Vmax)=max_Vmax;

Nrub_la=Nrub_la_mmol*mmoltogm;

%% Chlorophyll
Eta = Rnchl;
K_chl_lhc=NlhcF;

Phi = K_chl_lhc * K_lhc_leaf;

% Phi = K_lhc_leaf;

Chl = (10^3)* sum(Nleaf_mmol) * (Phi/Eta);
Chl_opt=Chl;

end

