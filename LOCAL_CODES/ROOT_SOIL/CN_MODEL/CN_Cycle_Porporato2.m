function [VARIABLES] = ...
    CN_Cycle_Porporato2(rootfr,PARAMS, SWITCHES, VARIABLES, FORCING, CONSTANTS, VERTSTRUC, STORAGE)
%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                              FUNCTION CODE                            %%
%%                     LONGWAVE ATTENUATION CALCULATION                  %%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   See Poporato et al., 2003; Quijano et al., 2013                       %
%     for equation # references                                           %
%-------------------------------------------------------------------------%
%   Created by   : Darren Drewry                                          %
%   Modified by  : Juan Quijano                                           %
%                : Dongkook Woo                                           %
%   Date         : January 11, 2014                                       %
%-------------------------------------------------------------------------%
%                                                                         %
%   Calculate Dynamics of Soil Carbon and Nitrogen
%
%   INPUTS:
%       INITIAL CONDITIONS:
%       sm = soil moisture profile [m^3 / m^3]
%       porsl = porosity profile [m^3 / m^3]
%       Ts = soil temperature profile [C]
%       Ts_max = Maximum temperature [C]. Used to compute environmental factors
%       qq = water fluxes between layers [mm / s]
%       layeruptake = Uptake from each layer [mm/s]
%       layeruptake_All = Uptake from each layer by each species separate [mm/s]
%       rootfr  = Root distribution of fine roots.

%       ku_Amm =
%       ku_Nit =
%       ki_Amm = partitioning coefficient for ammonium immobilization [m^3 d / gC]
%       ki_Nit = partitioning coefficient for nitrate immobilization [m^3 d / gC]

%-------------------------------------------------------------------------%
%%

% DE REFERENCE BLOCKS
% water and temperature in the soil
nspecies = PARAMS.CanStruc.nspecies; % number of species
if SWITCHES.CN.Bioturbation
    
    %     VARIABLES.SOIL.volliqli=VARIABLES.SOIL.volliq(1);
    
    Theta = [VARIABLES.SOIL.volliqli ; VARIABLES.SOIL.volliq];
    Ts = [VARIABLES.SOIL.Tli ; VARIABLES.SOIL.Ts];
    qq = [VARIABLES.SOIL.qinflL ; VARIABLES.SOIL.qlayer*ones(1,nspecies)];
    layeruptake_all_Sp = [zeros(1,nspecies) ; VARIABLES.SOIL.layeruptake_all];
    % Dongkook Woo - Edit
    %smp = [VARIABLES.SOIL.psil ; VARIABLES.SOIL.smp];
    smp = [VARIABLES.SOIL.psili ; VARIABLES.SOIL.smp];
    % Dongkook Woo - Edit End
    nl_soil=PARAMS.nl_soil + 1;%       nl_soil = # soil layers
    rootfr=[zeros(1,nspecies) ; rootfr];
    dz_mm = [VARIABLES.SOIL.litterthickness*1000 ; VERTSTRUC.dzsmm*ones(1,nspecies)];%       dz_mm = grid depth [mm]
else
    Theta = VARIABLES.SOIL.volliq;
    Ts = VARIABLES.SOIL.Ts;
    qq = VARIABLES.SOIL.qlayer*ones(1,nspecies);
    layeruptake_all_Sp = VARIABLES.SOIL.layeruptake_all;
    smp = VARIABLES.SOIL.smp;
    nl_soil=PARAMS.nl_soil ;%       nl_soil = # soil layers
    rootfr=rootfr;
    dz_mm = VERTSTRUC.dzsmm*ones(1,nspecies);%       dz_mm = grid depth [mm]
end

a_Amm = PARAMS.CN.a_Amm;%       a_Amm = fraction of dissolved ammonium [-]
a_Nit = PARAMS.CN.a_Nit;%       a_Nit = fraction of dissolved nitrate [-]
ki_Amm = PARAMS.CN.ki_Amm; % constant that determines the partition between NO3 and NH4 for immobilization
ki_Nit = PARAMS.CN.ki_Nit; % constant that determines the partition between NO3 and NH4 for immobilization
rr = PARAMS.CN.rr;%       rr = fraction of decomposed organic carbon that goes to respiration [-]

VARIABLES.Amm(isinf(VARIABLES.Amm)|isnan(VARIABLES.Amm)) = 0;
VARIABLES.Ammdd(isinf(VARIABLES.Ammdd)|isnan(VARIABLES.Ammdd)) = 0;
VARIABLES.Nit(isinf(VARIABLES.Nit)|isnan(VARIABLES.Nit)) = 0;
VARIABLES.Nitdd(isinf(VARIABLES.Nitdd)|isnan(VARIABLES.Nitdd)) = 0;



TR_can_store = VARIABLES.CANOPY.TR_can_all;
dtime = CONSTANTS.dtime;%       dtime = Time step [1800 s]
%*************************************************************************
%                       PREALLOCATE VECTORS
PHI = zeros(nl_soil,1);
phi = ones(nl_soil,1);
%*************************************************************************
%                       CHANGE OF UNITS AND DE-REFERENCE BLOCK

qq=qq*86400/1000;  % Change units from [mm/s] to [m/d]

%layeruptake=layeruptake*86400/1000;  % Change units from [mm/s] to [m/d]
layeruptake_all_Sp=layeruptake_all_Sp*86400/1000;  % Change units from [mm/s] to [m/d]

% check that evaporation fluxes are not included
if SWITCHES.CN.Bioturbation
    qq(1) = max(qq(1),0);
    qq(2) = max(qq(2),0);
end
if qq(2) > 0
    %qq(2) = 0;
    stophere=0;
end
%%qq(1)=0;
%qq(2)=0;
%*************************************************************************
%     Applying N deposition
%*************************************************************************

CN_envfactor ();

for si=1:nspecies
    
    sm=Theta(:,si);
    
    layeruptake_all=layeruptake_all_Sp(:,si);
    
    Cl = VARIABLES.Cl(:,si);%       Cl = carbon concentration in litter pool [gC / m^3]
    Ch = VARIABLES.Ch(:,si);%       Ch = carbon concentration in humus pool [gC / m^3]
    Cb = VARIABLES.Cb(:,si);%       Cb = carbon concentration in biomass pool [gC / m^3]
    Nl = VARIABLES.Nl(:,si);%       Nl = nitrogen concentration in litter pool [gN / m^3]
    
    Amm = VARIABLES.Amm(:,si);%       Amm = ammonium concentration in soil [gN / m^3]
    Ammdd = VARIABLES.Ammdd(:,si);%       Amm = ammonium concentration in zone of disturbance [gN / m^3]
    Nit = VARIABLES.Nit(:,si);%       Nit = nitrate concentration in soil [gN / m^3]
    Nitdd = VARIABLES.Nitdd(:,si);%       Nit = nitrate concentration in zone of disturbance [gN / m^3]
    
    
    
    if SWITCHES.CN.AtmosDeposition
        [VARIABLES] = CN_deposition (VARIABLES,PARAMS,CONSTANTS,VERTSTRUC, si);
        Amm = VARIABLES.Amm(:,si);
        Nit = VARIABLES.Nit(:,si);
        Ammdd = VARIABLES.Ammdd(:,si);
        Nitdd = VARIABLES.Nitdd(:,si);
    end
    
    
    %*************************************************************************
    %     Applying N Fertilizer
    %*************************************************************************
    if SWITCHES.CN.Fertilizer(si)
        if FORCING.doy==PARAMS.CN.N_Fert_DOY
            disp('Fertilizer applied')
        end
        if floor(FORCING.doy)==PARAMS.CN.N_Fert_DOY
            [VARIABLES] = CN_fertilizer (VARIABLES,PARAMS,CONSTANTS,VERTSTRUC,si);
            Amm = VARIABLES.Amm(:,si);
            Nit = VARIABLES.Nit(:,si);
            Ammdd = VARIABLES.Ammdd(:,si);
            Nitdd = VARIABLES.Nitdd(:,si);
        end
    end
    
    %*************************************************************************
    %     COMPUTE ENVIRONMENTAL FACTORS
    %*************************************************************************
    
    %*************************************************************************
    %   BIOTURBATION (take out)
    %*************************************************************************
    if SWITCHES.CN.Bioturbation
        [Cl, VARIABLES] = CN_bioturbation (PARAMS, VARIABLES, CONSTANTS,...
            FORCING, VERTSTRUC, SWITCHES, fTd, si);
    end
    
    
    
    
    %*************************************************************************
    %     COMPUTE ADITTION OF LITTER. (LITTER, HUMUS, BIOMASS)
    %*************************************************************************
    [CNa, ADD, CNo, OUT, ADD_bio, ADD_ex, ADD_net] = CN_addlitter (FORCING, PARAMS,...
        VERTSTRUC, VARIABLES, CONSTANTS, SWITCHES, rootfr, si);
    
    %*************************************************************************
    %     COMPUTE RATES OF CHANGES IN POOLS CONCENTRATION. (LITTER, HUMUS, BIOMASS)
    %*************************************************************************
    [dCl_dt, dNl_dt, dCh_dt, dNh_dt, dCb_dt, dNb_dt, DECl, DECh, VARIABLES] = ...
        CN_dynamics(VARIABLES, PARAMS, SWITCHES, fSd, fTd, phi, PHI, ADD, CNa, OUT, CNo, si);
    
    
    %*************************************************************************
    %     COMPUTE NET MINERALIZATION OR NET IMMOBILIZATION
    %*************************************************************************
    % Preallocate vectors
    [phi, PHI, MIN_net, IMM_net, MIN_gross, IMM_gross, Nreg, DECl, DECh] = ...
        CN_computephi (VARIABLES, PARAMS, SWITCHES, fSd, fTd, phi, ADD, CNa, si);
    
    
    [dCl_dt, dNl_dt, dCh_dt, dNh_dt, dCb_dt, dNb_dt, DECl, DECh, VARIABLES] =...
        CN_dynamics (VARIABLES, PARAMS, SWITCHES, fSd, fTd, phi, PHI, ADD, CNa, OUT, CNo, si);
    
    
    
    % PARTITION IMMOBILIZATION BETWEEN NO3- AND NH4+
    indn = Amm ~= 0;
    IMM_Amm(indn) = (ki_Amm*Amm(indn)./ (ki_Amm*Amm(indn) + ki_Nit*Nit(indn))).*IMM_net(indn);
    IMM_Amm(~indn) = 0;
    IMM_Amm = IMM_Amm(:);
    indn = Nit ~= 0;
    IMM_Nit(indn) = (ki_Nit*Nit(indn)./ (ki_Amm*Amm(indn) + ki_Nit*Nit(indn))).*IMM_net(indn);
    IMM_Nit(~indn) = 0;
    IMM_Nit = IMM_Nit(:);
    
    IMM_net_all(:,si)=IMM_net;
    
    DECl_all(:,si)=DECl;
    
    DECh_all(:,si)=DECh;
    
    MIN_net_all(:,si)=MIN_net;
    IMM_Amm_all(:,si)=IMM_Amm;
    IMM_Nit_all(:,si)=IMM_Nit;
    
    
    dCl_dt_all(:,si)=dCl_dt;
    dCh_dt_all(:,si)=dCh_dt;
    dCb_dt_all(:,si)=dCb_dt;
%     dNh_dt_all(:,si)=dNh_dt;
%     dNb_dt_all(:,si)=dNb_dt;
    dNl_dt_all(:,si)=dNl_dt;
    
    MIN_gross_all(:,si)=MIN_gross;
    IMM_gross_all(:,si)=IMM_gross;
    
    Nreg_all(:,si)=Nreg;
    
    ADD_all(:,si)=ADD;
    
    ADD_net_all(:,si)=ADD_net;
    
        
end
    
Amm_all = VARIABLES.Amm;
Nit_all = VARIABLES.Nit;
Ammdd_all = VARIABLES.Ammdd;
Nitdd_all = VARIABLES.Nitdd;


Cl_all = VARIABLES.Cl;
Ch_all = VARIABLES.Ch;
Cb_all = VARIABLES.Cb;
Nl_all = VARIABLES.Nl;


if SWITCHES.CN.NupRootBiomass == 1
    %*************************************************************************
    %     COMPUTE NITRATE AND AMMONIUM LEACHING FLUX [gN / m^3 / d]
    %*************************************************************************
    [rootarq] = CN_rootprop (FORCING, PARAMS, VARIABLES, VERTSTRUC, SWITCHES);
    RVID = sum(rootarq.RVID,2);
    
    % update concentrations in Nit and Amm by root death or root growth
    if VARIABLES.timestep > 1          % timestep = Current time step
        CN_updatecon(VERTSTRUC, SWITCHES, VARIABLES, PARAMS, rootarq.RVID, Nit_all, Nitdd_all);
    end
    
    [VARIABLES] = CN_nit_transport(PARAMS, SWITCHES, VARIABLES, VERTSTRUC, CONSTANTS,...
        FORCING, rootarq, Nit_all, Nitdd_all, a_Nit, qq, Theta, layeruptake_all_Sp);
    
    [VARIABLES] = CN_amm_transport(PARAMS, SWITCHES, VARIABLES, VERTSTRUC, CONSTANTS,...
        FORCING, rootarq, Amm_all, Ammdd_all, a_Amm, qq, Theta, layeruptake_all_Sp);
    
    % store leaching in units of [gr / m^3 / d]
    LCH_N_m2 = VARIABLES.LCH_nit_m2 + VARIABLES.LCH_amm_m2;  %[gr/m2/d]
    % Dongkook - Edit
    LCH_nit_m2 = VARIABLES.ADV_nit_m2+VARIABLES.DIFF_nit_m2;
    LCH_amm_m2 = VARIABLES.ADV_amm_m2+VARIABLES.DIFF_amm_m2;
    LCH_nit = LCH_nit_m2./(dz_mm/1000);
    LCH_amm = LCH_amm_m2./(dz_mm/1000);
    % Dongkook - Edit End
    % Compute total Nitrogen Uptake [gr/m2/d].
    
    UP_N_m2 = sum(sum(VARIABLES.UP_nit_m2) + sum(VARIABLES.UP_amm_m2));
    
    UP_N_m2_sp = sum(VARIABLES.UP_nit_m2) + sum(VARIABLES.UP_amm_m2);   %[gr/m^2/d]
    
elseif SWITCHES.CN.NupRootBiomass == 0
    %*************************************************************************
    %     COMPUTE NITRATE AND AMMONIUM LEACHING FLUX [gN / m^3 / d]
    %*********************************************R****************************
    [LCH_nit_m2, TLCH_nit_m2] = CN_leach(PARAMS, SWITCHES, Nit_all, a_Nit, qq, Theta);
    [LCH_amm_m2, TLCH_amm_m2] = CN_leach(PARAMS, SWITCHES, Amm_all, a_Amm, qq, Theta);
    % Compute leaching in units of [gr / m^3 / d]
    LCH_nit = LCH_nit_m2./(dz_mm/1000);   % [gr/m^3/d]
    LCH_amm = LCH_amm_m2./(dz_mm/1000);   % [gr/m^3/d]
    
    LCH_N_m2 = TLCH_nit_m2 + TLCH_amm_m2;  %[gr/m2/d] % Out of system
    %LCH_N_m2 = sum((LCH_nit+LCH_amm).*(dz_mm/1000));
    %LCH_N_m2 = sum(LCH_amm_m2+LCH_nit_m2);
    %*************************************************************************
    %     COMPUTE NITRATE AND AMMONIUM WATER UPTAKE FLUX [gN / m^2 / d]
    %*************************************************************************
    %   [SUP_amm_m2, SUP_nit_m2, SUP_amm_all_m2, SUP_nit_all_m2, SUP_amm_active_m2, SUP_nit_active_m2, SUP_amm_passive_m2, SUP_nit_passive_m2,...
    %       Amm_DEM, Nit_DEM, Fix_Nit_m2, Fix_Amm_m2, DEM_fraction_l_s_g_b, VARIABLES, UPr_amm_m2, UPr_nit_m2] = CN_nuptake (VARIABLES, PARAMS, SWITCHES, VERTSTRUC,...
    %       a_Nit, a_Amm, sm, layeruptake_all, CARBONS, CNratio_leaf, CNratio_stem, CNratio_grain, CNratio_root, rootfr, choose_crop,...
    %       scenarios, totlail_layer_day_point, Add_Carbon, how_many_year,...
    %       detect_first_day_of_year, day, iscorn,...
    %       DEM_fraction_l_s_g_b_rate,Check_year,Tot_UP_Amm_m2_print,Tot_UP_Nit_m2_print);
    
    [VARIABLES, SUP_amm_pas_all_m2, SUP_nit_pas_all_m2, SUP_amm_act_all_m2, SUP_nit_act_all_m2, ...
        FUP_amm_all_m2, FUP_nit_all_m2, RUP_amm_all_m2, RUP_nit_all_m2] = ...
        CN_nuptake2 (VARIABLES, PARAMS, SWITCHES, VERTSTRUC, FORCING,...
        STORAGE, CONSTANTS, a_Nit, a_Amm, Theta, layeruptake_all_Sp, rootfr);
    
    %     % Compute uptake in units of [gr / m^3 / d]
    %     UP_nit = UP_nit_m2./(dz_mm/1000);   % [gr/m^3/d]
    %     UP_amm = UP_amm_m2./(dz_mm/1000);   % [gr/m^3/d]
    %     UP_nit_all = UP_nit_all_m2./repmat((dz_mm/1000),1,nspecies); % [gr/m^3/d]
    %     UP_amm_all = UP_amm_all_m2./repmat((dz_mm/1000),1,nspecies); % [gr/m^3/d]
    TotUP_nit_all = (SUP_nit_pas_all_m2+SUP_nit_act_all_m2+FUP_nit_all_m2+RUP_nit_all_m2)./(dz_mm/1000); % [gr/m^3/d]
    TotUP_amm_all = (SUP_amm_pas_all_m2+SUP_amm_act_all_m2+FUP_amm_all_m2+RUP_amm_all_m2)./(dz_mm/1000); % [gr/m^3/d]
    
    SUP_nit_all = (SUP_nit_pas_all_m2+SUP_nit_act_all_m2)./(dz_mm/1000); % [gr/m^3/d]
    SUP_amm_all = (SUP_amm_pas_all_m2+SUP_amm_act_all_m2)./(dz_mm/1000); % [gr/m^3/d]
    SUP_nit = sum(SUP_nit_all,2); % [gr/m^3/d]
    SUP_amm = sum(SUP_amm_all,2); % [gr/m^3/d]
    
    FUP_nit_all = (FUP_nit_all_m2)./(dz_mm/1000); % [gr/m^3/d]
    FUP_amm_all = (FUP_amm_all_m2)./(dz_mm/1000); % [gr/m^3/d]
    
    RUP_nit_all = (RUP_nit_all_m2)./(dz_mm/1000); % [gr/m^3/d]
    RUP_amm_all = (RUP_amm_all_m2)./(dz_mm/1000); % [gr/m^3/d]
    
    % Compute total Nitrogen Uptake from soil [gr/m2/d].
    %UP_N_m2 = sum(UP_amm_m2) + sum(UP_nit_m2);   %[gr/m^2/d]
    UP_N_m2 = sum(sum(SUP_nit_pas_all_m2+SUP_nit_act_all_m2+SUP_amm_pas_all_m2+SUP_amm_act_all_m2)); %[gr/m^2/d]

    UP_N_m2_sp= sum(SUP_nit_pas_all_m2+SUP_nit_act_all_m2+SUP_amm_pas_all_m2+SUP_amm_act_all_m2); %[gr/m^2/d]

end

 



%% ***************************************


%*************************************************************************
%       CHANGES IN MINERAL CONCENTRATIONS
%*************************************************************************
%    [Nitrif, Denitrif, Volat] = CN_nfluxes (VARIABLES, PARAMS, VERTSTRUC, SWITCHES, fSn, fTn, DECl);
[Nitrif, Denitrif, Volat, Total_N20, N2O_nit] = CN_nfluxes (VARIABLES,...
    PARAMS, VERTSTRUC, SWITCHES, fSn, fTn, DECl_all, DECh_all);
%*************************************************************************
%      STATE UPDATES
%*************************************************************************
% nitrate and ammonium in the soil
if SWITCHES.CN.NupRootBiomass == 1
    dAmm_dt = MIN_net_all - IMM_Amm_all - Nitrif - Volat + VARIABLES.deltaamm;  %[gr N / m^3 / d]
    dNit_dt = Nitrif - IMM_Nit_all - Denitrif + VARIABLES.deltanit;%[gr N / m^3 / d]
elseif SWITCHES.CN.NupRootBiomass == 0
    dAmm_dt = MIN_net_all - IMM_Amm_all - Nitrif - Volat - SUP_amm_all - LCH_amm; %  [gr N / m^3 / d]
    dNit_dt = Nitrif - IMM_Nit_all - Denitrif - SUP_nit_all - LCH_nit; %[gr N / m^3 / d]
end
if SWITCHES.CN.Denitrification == 1
    dNit_dt = dNit_dt-N2O_nit;
else
    dNit_dt = dNit_dt;
end

% nitrate and ammonium in the zs
if SWITCHES.CN.NupRootBiomass == 1
    dAmmdd_dt = VARIABLES.deltaammdd;
    dNitdd_dt = VARIABLES.deltanitdd;
end

Cl_new = Cl_all + dCl_dt_all*dtime/86400;          %[gC/m3]
if SWITCHES.CN_type
    Ch_new = Ch_all + dCh_dt_all*dtime/86400;      %[gC/m3]
    % Dongkook - Edit
    %        Nh=VARIABLES.Nh;
    %        Nh_new = Nh + dNh_dt*dtime/86400;
    % Dongkook - Edit End
else
    Ch_new = nan*(zeros(nl_soil,1));       %[gC/m3]
end

%     dAmm_dt=0.1;
%     dNl_dt=0.1;

Cb_new = Cb_all + dCb_dt_all*dtime/86400;          %[gC/m3]
Nl_new = Nl_all + dNl_dt_all*dtime/86400;          %[gC/m3]
Amm_new = Amm_all + dAmm_dt*dtime/86400;       %[gC/m3]
indneg = Amm_new < 0;

if sum(indneg)==1
    stop =4;
end
Amm_new(indneg) = 0;                       %[gC/m3]

Nit_new = Nit_all + dNit_dt*dtime/86400;       %[gC/m3]
indneg = Nit_new < 0;
if sum(indneg)==1
    stop =4;
end
Nit_new(indneg) = 0;
CNl_new = Cl_new./Nl_new;

indcb = Cb_new<1;
Cb_new(indcb) = 1;
%*************************************************************************
% STORE IN STRUCTURE
%*************************************************************************
CN_massbalance ();


%*************************************************************************
% STORE IN STRUCTURE
%*************************************************************************

VARIABLES.Cl=real(Cl_new);
VARIABLES.Ch=real(Ch_new);
VARIABLES.Cb=real(Cb_new);
VARIABLES.Nl=real(Nl_new);
VARIABLES.Amm=real(Amm_new);
VARIABLES.Nit=real(Nit_new);
VARIABLES.CNl=real(CNl_new);
VARIABLES.dCl_dt = real(dCl_dt_all);
VARIABLES.dCh_dt = real(dCh_dt_all);
VARIABLES.dCb_dt = real(dCb_dt_all);

% CARBON MASS BALANCE
VARIABLES.Cinput = real(Cinput); %[grC/m2/d]  Input from Litter and Bioturbation
VARIABLES.CBoutput = real(CBoutput);  %[grC/m2/d]  Output from Bioturbation
VARIABLES.DECout = real(DECout);   %[grC/m2/d] Respiration
VARIABLES.dCb_m2 = real(dCb_m2);   %[grC/m2/d] change in Microbial Biomass
VARIABLES.dCl_m2 = real(dCl_m2);   %[grC/m2/d] change in organic matter
VARIABLES.mberrorC = real(mberrorC); %[grC/m2/d]

% NITROGEN MASS BALANCE
VARIABLES.UP_N_m2=real(UP_N_m2);     %[grC/m2/d]

VARIABLES.UP_N_m2_sp=real(UP_N_m2_sp);

if SWITCHES.CN.NupRootBiomass == 0
    VARIABLES.UP_nit_m3 = sum(TotUP_nit_all,2);                    %[g/m^3/d]
    VARIABLES.UP_nit_m2 = sum(TotUP_nit_all,2).*(dz_mm/1000);                    %[g/m^2/d]
    VARIABLES.UP_nit_all_m2 = (SUP_nit_pas_all_m2+SUP_nit_act_all_m2+FUP_nit_all_m2+RUP_nit_all_m2);
    VARIABLES.UP_amm_m3 = sum(TotUP_amm_all,2);                     %[g/m^3/d]
    VARIABLES.UP_amm_m2 = sum(TotUP_amm_all,2).*(dz_mm/1000);                   %[g/m^2/d]
    VARIABLES.UP_amm_all_m2 = (SUP_amm_pas_all_m2+SUP_amm_act_all_m2+FUP_amm_all_m2+RUP_amm_all_m2);
end
if SWITCHES.CN.NupRootBiomass == 1
    VARIABLES.dAmmdd_m2 = real(dAmmdd_m2); %[grC/m2/d]
    VARIABLES.dNitdd_m2 = real(dNitdd_m2); %[grC/m2/d]
end

VARIABLES.MIN_net=real(MIN_net_all);
VARIABLES.MIN_gross=real(MIN_gross_all);
VARIABLES.LCH_N_m2=real(LCH_N_m2);  %[grC/m2/d]
VARIABLES.LCH_nit_m2=real(LCH_nit_m2);
VARIABLES.LCH_nit_m3=real(LCH_nit);
VARIABLES.LCH_amm_m2=real(LCH_amm_m2);
VARIABLES.LCH_amm_m3=real(LCH_amm);

VARIABLES.Ninput = real(Ninput);  %[grC/m2/d]
VARIABLES.NBoutput = real(NBoutput); %[grC/m2/d]
VARIABLES.dNb_m2 = real(dNb_m2); %[grC/m2/d]
VARIABLES.dNl_m2 = real(dNl_m2); %[grC/m2/d]
VARIABLES.dAmm_m2 = real(dAmm_m2); %[grC/m2/d]
VARIABLES.dNit_m2 = real(dNit_m2); %[grC/m2/d]
VARIABLES.IMM_net=real(IMM_net_all);
VARIABLES.IMM_gross=real(IMM_gross_all);
VARIABLES.Nreg=real(Nreg_all);
VARIABLES.DECl=real(DECl_all);
VARIABLES.PHI=real(PHI);
VARIABLES.phi=real(phi);
VARIABLES.fSd=real(fSd);
VARIABLES.fTd=real(fTd);
VARIABLES.mberrorN=real(mberrorN);
VARIABLES.mberrorC=real(mberrorC);

VARIABLES.ADD = real(ADD_all);
VARIABLES.ADD_bio = real(ADD_bio);
VARIABLES.ADD_ex = real(ADD_ex);
VARIABLES.OUT = real(OUT);
VARIABLES.ADD_net = real(ADD_net_all);

if SWITCHES.CN.NupRootBiomass == 1
    VARIABLES.RMI = rootarq.RMI;
    VARIABLES.RLI = rootarq.RLI;
    VARIABLES.RAI = rootarq.RAI;
    VARIABLES.RVI = rootarq.RVI;
    VARIABLES.RAID = rootarq.RAID;
    VARIABLES.RVID = rootarq.RVID;
end
% Dongkook Woo - Edit
% Denitrification
VARIABLES.Denitrif=real(Denitrif);
VARIABLES.Total_N20=real(Total_N20);
VARIABLES.N2O_nit=real(N2O_nit);

VARIABLES.Nitrif = real(Nitrif);
VARIABLES.DECh = real(DECh_all);
% Dongkook Woo - Edit End
% % 
% if VARIABLES.timestep == 6621
%     PPPPPPPP;
% end

