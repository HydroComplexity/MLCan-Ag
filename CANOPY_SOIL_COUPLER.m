
nl_can = PARAMS.CanStruc.nl_can;
nl_soil=PARAMS.nl_soil;


for si=1:num_species
    file_name_Py   = sprintf('YrwdSpec_%sY%s.mat', num2str(si), num2str(Start_Y));
    
    load('.\Temps\temporary.mat', 'config_path');
    
    YrPath=fullfile(config_path, 'YearFr');
    
    fullpath_Py    =...
        fullfile(YrPath, file_name_Py);
    
    load (fullpath_Py,...
        'dat_LAD')
    
    znc = (dat_LAD(:,1))';                                               % Middle of each canopy grid [m]
    zhc = (dat_LAD(:,1)+(dat_LAD(2,1)-dat_LAD(1,1))/2)';               % Top of each canopy grid [m]
    %dzc = zhc(1)./2;
    dzc = zhc(1);
    
    VERTSTRUC.znc = znc;
    VERTSTRUC.zhc = zhc;
    VERTSTRUC.dzc = dzc;
    
    LADnorm_all(:,si)= dat_LAD(:,2);
    
    
    %% Soil Charecteristics
    
    load (fullpath_Py,...
        'dat_root')
    
    zns = dat_root(:,1);
    dzs=zeros(nl_soil,1);
    % Soil layer thicknesses
    dzs(1)  = 0.5*(zns(1)+zns(2));
    dzs(nl_soil)= zns(nl_soil)-zns(nl_soil-1);
    for j = 2:nl_soil-1
        dzs(j)= 0.5*(zns(j+1)-zns(j-1));
    end
    
    %     dzs=dzs';
    
    zhs=zeros(nl_soil,1);
    
    % Soil layer interface depths from the surface [m]
    zhs(nl_soil) = zns(nl_soil) + 0.5*dzs(nl_soil);
    for j = 1:nl_soil-1
        zhs(j)= 0.5*(zns(j)+zns(j+1));
    end
    %     zhs=zhs';
    
    znsmm = zns(:)*1000;      % [mm]
    dzsmm = dzs(:)*1000;      % [mm]
    zhsmm = zhs(:)*1000;      % [mm]
    
    rootfr(:,si) = dat_root(:,2);
    roottr(:,si) = dat_root(:,2);
    
    count_nl_root=dat_root(:,2);
    count_nl_root(1)=1;
    nl_root(si) = sum(count_nl_root ~= 0);
    
end

% INITIALIZE CANOPY STATES


% ASSIGN
VERTSTRUC.zns = zns;
VERTSTRUC.dzs = dzs;
VERTSTRUC.zhs = zhs;
VERTSTRUC.znsmm = znsmm;
VERTSTRUC.dzsmm = dzsmm;
VERTSTRUC.zhsmm = zhsmm;
VERTSTRUC.rootfr = rootfr;
VERTSTRUC.roottr = roottr;
VERTSTRUC.nl_root = nl_root;
VERTSTRUC.nl_soil = nl_soil;
PARAMS.Soil.nl_soil = nl_soil;

[VERTSTRUC] = SOIL_PROPERTIES(PARAMS, VERTSTRUC);
theta_dry = VERTSTRUC.theta_dry;
porsl = VERTSTRUC.porsl;
TK_dry = VERTSTRUC.TK_dry;
TK_sol = VERTSTRUC.TK_sol;
HC_sol = VERTSTRUC.HC_sol;

%%
% INITIALIZE SOIL STATES
% initialize soil moisture
% if sum(volliqinit>=porsl) > 0
%     volliqinit= porsl;
% end
% if sum(volliqinit<=theta_dry) > 0
%     volliqinit= theta_dry;
% end
for ss=1:nl_soil
    if volliqinit(ss)>=porsl(ss)
        volliqinit(ss)= porsl(ss);
    end
    if volliqinit(ss)<=theta_dry(ss)
        volliqinit(ss)= theta_dry(ss);
    end
end

VARIABLES.SOIL.volliq = volliqinit.*ones(1, num_species);
VARIABLES.SOIL.volliqli = volliqliinit*ones(1, num_species) ;   % Initial value of litter soil moisture

% initialize ice content & soil temperature
Ts = Tsinit*ones(1, num_species);
VARIABLES.SOIL.Tli = Ta_in(1)*ones(1, num_species);
VARIABLES.SOIL.Tsl = Tslint*ones(1, num_species);
VARIABLES.SOIL.Tlprev = Ta_in(1)*ones(1, num_species);
VARIABLES.SOIL.TKsoil=VERTSTRUC.TK_sol*ones(1, num_species);


% RUN CANOPY-ROOT-SOIL MODEL
% LOOP OVER EACH YEAR TO RE-INITIALIZE CANOPY/SOIL STATES FOR EACH YEAR


%% Rohit Start
%############ Crop Growth Model: Initialize the variables #############%
% Num_years = CROP_GROWTH.Num_years;

dtime=CONSTANTS.dtime;

NoDay_prev  = 0;

%% RohitN end

for yy = 1:length(Run_years)
    % compute the range of time steps in current year
    yy;
    
    VARIABLES.yy=yy;
    
    % TimeBar 1/3
    tff=yendinds(yy);
    t00=ybeginds(yy);
    hh = timebar(['Progress:', num2str(Run_years(yy)),'/', num2str(Run_years(end))],'MLCan Simulation');
    tic
    
    ybind = ybeginds(yy);
    yeind = yendinds(yy);
    
    disp (['Year: ', num2str(Run_years(yy)), '/', num2str(Run_years(end))]);
    
    %% Rohit Start
    %####################################%
    if(rem(Run_years(yy), 4)== 0)
        NoDay=366;
    else
        NoDay=365;
    end
    
    VARIABLES.NoDay=NoDay;
    
    
    % ALLOCATE STORAGE FOR MODELLED VARIABLES
    ALLOCATE_STORAGE;
    
    %     NoDay=yeind/24;
    
    VARIABLES.CANOPY.gsv_sun = 0.01*ones(nl_can,nspecies);
    VARIABLES.CANOPY.gsv_shade = 0.01*ones(nl_can,nspecies);
    VARIABLES.CANOPY.TR = zeros(length(znc),nspecies);
    VARIABLES.CANOPY.Sh2o_prof = zeros(length(znc),1);
    VARIABLES.CANOPY.Tl_prev_dt = Ta_in(1) * ones(nl_can,1);
    
    
    VARIABLES.SOIL.smp = VERTSTRUC.psi0 .* (VARIABLES.SOIL.volliq ./ VERTSTRUC.porsl).^(-VERTSTRUC.bsw);
    
    % Message box for soil moisture
    if sum(volliqinit>porsl) > 0
        msgbox({'Initial soil moisture is higher than saturated soil moisture!', 'Solution: Modify initilal soil moisture or Increase % of sand.'},'error');
        error('Initial soil moisture is higher than saturated soil moisture!')
    end
    if sum(volliqinit<theta_dry) > 0
        msgbox({'Initial soil moisture is lower than residual soil moisture!', 'Solution: Modify initilal soil moisture or Increase % of sand.'},'error');
        error('Initial soil moisture is lower than residual soil moisture!')
    end
    
    % Initialize snow moisture variables
    VARIABLES.SOIL.voltotsn = zeros(1, num_species);
    VARIABLES.SOIL.voltotli = zeros(1, num_species);
    VARIABLES.SOIL.voliceli = zeros(1, num_species);
    VARIABLES.SOIL.volliqsn = zeros(1, num_species);
    VARIABLES.SOIL.volicesn = zeros(1, num_species);
    
    VARIABLES.SOIL.zicesl = zeros(1, num_species);
    
    VARIABLES.SOIL.zicesl_prev= VARIABLES.SOIL.zicesl;
    
    VARIABLES.SOIL.wicesl = zeros(1, num_species);
    
    VARIABLES.SOIL.rhosn = 1000*ones(1, num_species);
    
    % Fixing the constant soil layer problem
    VARIABLES.SOIL.volice = zeros(nl_soil, num_species);
    VARIABLES.SOIL.snow_tcount = zeros(1, num_species);
    
    
    VARIABLES.CANOPY.Sh2o_can_prev = 0;
    
    
    % INITIALIZE ROOT POTENTIAL
    for si=1:nspecies
        VARIABLES.ROOT.rpp_wgt(:,si) =  VARIABLES.SOIL.smp(1);
    end
    
    VARIABLES.ROOT.rpp= VARIABLES.SOIL.smp;
    
    % PEDOTRANSFER FUNCTIONS
    if SWITCHES.Pedofunctions
        [VERTSTRUC] = PEDOSOIL_PROPERTIES(PARAMS, VERTSTRUC, VARIABLES);
        porsl = VERTSTRUC.porsl;                                         % POROSITY
        psi0 = VERTSTRUC.psi0;                                           % MINIMUM SOIL SUCTION = SOIL POTENTIAL AT SATURATION [mm]
        bsw = VERTSTRUC.bsw;                                             % B PARAMETER BROKS AND COREY SHAPE PARAMETER
        Ksat = VERTSTRUC.HKsat;                                          % HYDRAULIC CONDUCTIVITY AT SATURATION [mm / s]
        eff_poros = VERTSTRUC.eff_poros;
    end
    
    StepsT    = CONSTANTS.timestep;
    
    Ta_tot=Ta_crop(ybind:yeind);
    
    PPT_tot=PPT_crop(ybind:yeind);
    
    [Max_Ta, Min_Ta, daily_ppt] = Convert_to_Daily(Ta_tot, NoDay, PPT_tot, StepsT);
    
    Hours_Temp=hour(ybind:yeind);
    intdoy = floor(round((doy_crop+1).*10^10)./10^10);
    
    
    %     Initialize Model Parameters
    InitializePar;
    
    
    %###################################%
    %% Rohit End
    
    %% LOOP OVER EACH TIME PERIOD IN YEAR yy
    %for tt = ybind:yeind
    for tt = ybind:1:yeind
        
        tt;
        % TimeBar 2/3
        timebar(hh,(tt-t00)/(tff-t00))
        
        timestep = tt-ybind + 1;
        VARIABLES.timestep = timestep;
        
        [VERTSTRUC, VARIABLES, rootfr] = ROOT_RESPONSE_DRY(VARIABLES,...
            SWITCHES, VERTSTRUC, CONSTANTS, PARAMS, doy, smp_store);
        
        % FORCING CONDITIONS
        FORCING.doy = doy(tt);
        FORCING.Rg = Rg_in(tt);
        FORCING.Pa = Pa_in(tt);
        if PARAMS.LWcom == 1
            FORCING.LWdn = LWdn_in(tt);
        end
        FORCING.zen = ZEN_in(tt);
        FORCING.U = U_in(tt);
        FORCING.ppt = PPT_in(tt);    % [mm]
        FORCING.Ta = Ta_in(tt);
        FORCING.ea = ea_in(tt);
        FORCING.Ca = CO2base;
        FORCING.ELEV=ELEV;
        CROP_GROWTH.LAI_actual=LAI_in(tt);
        
        if (~SWITCHES.soilheat_on)
            VARIABLES.SOIL.Ts = (Ta_in(tt)-5)*ones(nl_soil, num_species);
        else
            VARIABLES.SOIL.Ts = Ts;
        end
        
        VARIABLES.SOIL.Ts = Ts;
        VARIABLES.SOIL.Tsurf=Ts(1);
        
        
        tdt = 24 * 60 / CONSTANTS.timestep;
        tqt = tt + (CROP_GROWTH.DOY_start-1) * tdt;
        
        
        
        % %         CANOPY STRUCTURE
        for si=1:1:nspecies
            
            if SWITCHES.plants(si)
                
                if SWITCHES.LT == 1
                    LAILT = LAI_in(timestep,si);
                    if SWITCHES.CGM ==1
                        VERTSTRUC.LAIzall(:,si) = CROP_GROWTH.LAIsim_can(si)*LADnorm_all(:,si);
                    else
                        VERTSTRUC.LAIzall(:,si) = LAILT*LADnorm_all(:,si);
                    end
                else
                    if SWITCHES.CGM ==1
                        VERTSTRUC.LAIzall(:,si) = CROP_GROWTH.LAIsim_can(si).*LADnorm_all(:,si);
                    else
                        VERTSTRUC.LAIzall(:,si) = LAI_in(tt,si)*LADnorm_all(:,si);
                    end
                end
                
            else
                VERTSTRUC.LAIzall(:,si) = zeros(nl_can,1);
            end
        end
        
        
        if SWITCHES.CGM ==1
            LADnorm = sum(VERTSTRUC.LAIzall,2)/sum(CROP_GROWTH.LAIsim_can(:));
        else
            LADnorm = sum(VERTSTRUC.LAIzall,2)/sum(LAI_in(tt,:));
        end
        %
        LADnorm(isnan(LADnorm)) = 0;
        
        %% Rohit end
        
        %%
        VERTSTRUC.LAIz = sum(VERTSTRUC.LAIzall,2);
        VERTSTRUC.LADz = VERTSTRUC.LAIz ./ dzc; % Total LAD distribution
        fLAIz =VERTSTRUC.LAIzall./(repmat(sum(VERTSTRUC.LAIzall,2),1,nspecies));
        fLAIz(isnan(fLAIz)) = 0; % Set to zero whenever there is not LAI at a given layer
        VERTSTRUC.fLAIz = fLAIz; % Fraction of LAI in each species at each relative height level
        
        % create vinds
        % 1. For the total canopy
        LADmax = (max(VERTSTRUC.LAIzall,[],2)); % Maximum LAD
        
        nvinds = find(LADmax<=0);
        vinds = find(LADmax>0);
        
        VERTSTRUC.vinds = vinds;
        VERTSTRUC.nvinds = nvinds;
        
        %         LAI_ref_vert=CROP_GROWTH.LAI_ref(tt)*LADnorm_all;
        
        % 2. For All the species
        for si=1:nspecies
            nvinds_all{si} = find(VERTSTRUC.LAIzall(:,si) <= 0);
            vinds_all{si} = find(VERTSTRUC.LAIzall(:,si) > 0);
        end
        
        VERTSTRUC.nvinds_all = nvinds_all;
        VERTSTRUC.vinds_all = vinds_all;
        
        % INITIALIZE CANOPY ENVIRONMENT
        VARIABLES.CANOPY.TAz = Ta_in(tt) * ones(nl_can,nspecies);
        VARIABLES.CANOPY.CAz = CO2base * ones(nl_can,nspecies);
        VARIABLES.CANOPY.EAz = ea_in(tt) * ones(nl_can,nspecies);
        VARIABLES.CANOPY.PAz = Pa_in(tt) * ones(nl_can,1);
        VARIABLES.CANOPY.Uz = U_in(tt) * ones(nl_can,1);
        
        VARIABLES.CANOPY.TR_sun = zeros(nl_can,nspecies);
        VARIABLES.CANOPY.TR_shade = zeros(nl_can,nspecies);
        
        % INITIALIZE CANOPY STATES
        VARIABLES.CANOPY.Tl_can_sun = VARIABLES.CANOPY.TAz;
        VARIABLES.CANOPY.Tl_can_shade = VARIABLES.CANOPY.TAz;
        VARIABLES.CANOPY.Tl_sun = repmat(VARIABLES.CANOPY.TAz,1,nspecies);
        VARIABLES.CANOPY.Tl_shade = repmat(VARIABLES.CANOPY.TAz,1,nspecies);
        VARIABLES.CANOPY.Ci_sun = repmat(0.7 * VARIABLES.CANOPY.CAz,1,nspecies);
        VARIABLES.CANOPY.Ci_shade = repmat(0.7 * VARIABLES.CANOPY.CAz,1,nspecies);
        
        
        %% CANOPY MODEL SOLUTION
        [An_can, Ph_can, LE_can, H_can, dHcan, Rnrad_can, TR_can, ...
            Fc_soil, LE_soil, H_soil, Rnrad_soil, G, Tsurf, remainsoil,remaincan,remaineco,...
            Rnrad_sun, Rnrad_shade, Rnrad_eco, ...
            An_sun, An_shade, LE_sun, LE_shade, H_sun, H_shade, TR_sun, TR_shade, ...
            Tl_sun, Tl_shade, Tl_sun_up, Tl_shade_up, Tl_sun_lo, Tl_shade_lo, psil_sun, psil_shade, gsv_sun, gsv_shade, fsvg_sun, fsvm_sun,...
            fsvg_shade, fsvm_shade, Ci_sun, Ci_shade, CAz, TAz, EAz, Uz, gbv_sun, gbh_sun, gbv_shade, gbh_shade, ...
            LAI_sun, LAI_shade, fsun, fshade, ...
            Ph_limit_sun, Jc_C3_sun, Jj_C3_sun, Js_C3_sun, Jc_C4_sun, Jj_C4_sun, Js_C4_sun, ...
            Ph_limit_shade, Jc_C3_shade, Jj_C3_shade, Js_C3_shade, Jc_C4_shade, Jj_C4_shade, Js_C4_shade, ...
            PARabs_sun, PARabs_shade, NIRabs_sun, NIRabs_shade, SWout, ...
            LWabs_can, LWemit_soil, LWemit_can, LWemit_sun, LWemit_shade, LWout, LWoutM, RH_soil, fdiff, ...
            Sh2o_prof, Sh2o_can, ppt_ground, Ch2o_prof, Ch2o_can, Evap_prof, Evap_can, ...
            dryfrac, wetfrac, Vz, VARIABLES, FORCING,...
            SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out, SWsoildir_in, SWsoildir_out,...
            SWsoildif_in, SWsoildif_out, LWabs_canM, LWabs_soilM, LSshaCON, LSsunCON] = ...
            CANOPY_MODEL(SWITCHES, VERTSTRUC, FORCING, PARAMS, VARIABLES, CONSTANTS);
        
        An_can(An_can<-40)=-40;
        Ph_can(Ph_can<-40)=-40;
        
        %         VARIABLES.CANOPY.Ph_can_all
        
        CROP_GROWTH.ph_can = VARIABLES.CANOPY.Ph_can_all;
        
        
        %%         % Rohit Start
        
        %########## Using Crop Growth Model##########
        if SWITCHES.CGM ==1
            %-------------GDD and Onset selection (SM based)----------%
            
            for si=1:num_species
                
                %                 if SWITCHES.plants ==1
                
                TempBase    = PARAMS.CGM.GDDTBASE(si);
                TempCut     = PARAMS.CGM.GDDTCUT(si);
                
                [doy1, doy2, doy3, doy4, doy5, plantingDate, harvestingDate,...
                    CROP_GROWTH, VARIABLES] = Growing_Degree_Days(CROP_GROWTH,...
                    VARIABLES, CONSTANTS, PARAMS, SWITCHES, Max_Ta, Min_Ta, NoDay,...
                    TempCut, TempBase, yy, si, tt);
                
                CROP_GROWTH.plantingDate(si) = plantingDate;
                CROP_GROWTH.harvestingDate(si) = harvestingDate;
                
                CROP_GROWTH.doy1(si)=doy1;
                CROP_GROWTH.doy2(si)=doy2;
                CROP_GROWTH.doy3(si)=doy3;
                CROP_GROWTH.doy4(si)=doy4;
                CROP_GROWTH.doy5(si)=doy5;
                
                
                Tl_mean_cgm = VARIABLES.CANOPY.Tl_mean(si);
                ph_can_cgm= CROP_GROWTH.ph_can(si);
                
                %------------- Biomass Estimation -----------------%
                [CROP_GROWTH] = Biomass_Estimation( CROP_GROWTH, ...
                    CONSTANTS, PARAMS, ph_can_cgm, Tl_mean_cgm, doy1, ...
                    doy2, doy3, doy4, doy5, si, tt);
                
                
                DR = CROP_GROWTH.DeathRate(si)-VARIABLES.CANOPY.DeathRate_prev(si);
                RDR = CROP_GROWTH.RootDeathRate(si)-VARIABLES.CANOPY.RootDeathRate_prev(si);
                
                DR(DR<0)=0;
                RDR(RDR<0)=0;
                
                VARIABLES.CANOPY.DeathRate(si)=DR;
                
                VARIABLES.CANOPY.RootDeathRate(si)=RDR;
                
                %                 end
                
            end
            
            VARIABLES.CANOPY.DeathRate_prev = CROP_GROWTH.DeathRate;
            VARIABLES.CANOPY.RootDeathRate_prev = CROP_GROWTH.RootDeathRate;
            
        else
            
            VARIABLES.CANOPY.DeathRate = zeros(1,num_species);
            VARIABLES.CANOPY.RootDeathRate = zeros(1,num_species);
            
        end
        
        
        for si=1:num_species
            if tqt>=CROP_GROWTH.plantingDate(si)*tdt && tqt<=CROP_GROWTH.harvestingDate(si)*tdt
                
                SWITCHES.plants(si)=1;
                
            else
                
                SWITCHES.plants(si)=0;
                
            end
        end
               
        
        %% Rohit Start
        %########## Irrigation selection ##########%
        
        Irrigation_Selection;
        
        %##############################%
        
        %###########################%
        %% Rohit end
        
        
        %%
        % SOLUTION OF SNOW-LITTER PACK DYNAMICS
        
        for si=1:num_species
            if SWITCHES.litter(si)
                [VARIABLES] = FLUXES_WATER_SOIL_LITTER (PARAMS, VARIABLES, CONSTANTS, FORCING, SWITCHES, si);
            else
                [VARIABLES] = FLUXES_WATER_SOIL (PARAMS, VARIABLES, CONSTANTS, FORCING, SWITCHES, si);
            end
        end
        
        % assign
        %         qinfl = VARIABLES.SOIL.qinfl;
        %         qinflL = VARIABLES.SOIL.qinfl;
        %         net_qinflL = VARIABLES.SOIL.net_qinflL;
        %         drainlitter = VARIABLES.SOIL.drainlitter;
        %         volliqli = VARIABLES.SOIL.volliqli;
        
        
        % Implicit Solution
        if (SWITCHES.ns)
            if tt == 321
                stop = 1;
            end
            [rpp,rpp_wgt,krad,kax,dwat,smp,bk,hk, ...
                qlayer,layeruptake,layeruptake_all,mberrormm,type,...
                hor_drainage, hor_drainage_lay,flux_Ss]=ROOTSOIL(SWITCHES,...
                VERTSTRUC, PARAMS, VARIABLES, CONSTANTS, nspecies);
            
            VARIABLES.SOIL.flux_Ss =flux_Ss;
            VARIABLES.ROOT.rpp = rpp;
            VARIABLES.ROOT.rpp_wgt = rpp_wgt;
            VARIABLES.ROOT.krad = krad;
            VARIABLES.SOIL.type = type;
            VARIABLES.ROOT.kax= kax;
        else
            
            for si=1:num_species
                
                if (SWITCHES.HR_on)
                    [rpp, rpp_wgt, krad, kax] = ROOTS_HR(SWITCHES, VERTSTRUC, PARAMS, VARIABLES, si);
                else
                    [rpp, rpp_wgt, krad, kax] = ROOTS_NOHR(SWITCHES, VERTSTRUC, PARAMS, VARIABLES, si);
                end
                
                VARIABLES.ROOT.rpp(:,si) = rpp;
                VARIABLES.ROOT.rpp_wgt(si) = rpp_wgt;
                VARIABLES.ROOT.krad(:, si) = krad;
                VARIABLES.ROOT.kax(:, si)= kax;
                
                % Soil Moisture Solution
                [dwat, smp, hk, smp_wgt, thsatfrac_wgt, qlayer] = ...
                    SOILMOISTURE(SWITCHES, VERTSTRUC, PARAMS, VARIABLES, CONSTANTS, si);
                
            end
            
            VARIABLES.SOIL.type=[0,0,0];
            VARIABLES.SOIL.flux_Ss=zeros(nl_soil,1);
            
            layeruptake = (smp - VARIABLES.ROOT.rpp).*VARIABLES.ROOT.krad;
            layeruptake_all = layeruptake;
            hor_drainage = 0;
            hor_drainage_lay = zeros(length(smp),1);
            mberrormm = nan;
        end
        %
        % Update Volumetric Liquid Content
        volliq = VARIABLES.SOIL.volliq;
        
        for si=1:num_species
            
            volliq(:,si) = volliq(:,si) + dwat(:,si);
            volliq(:,si) = max(theta_dry, volliq(:,si));
            volliq(:,si) = min(VERTSTRUC.eff_poros, volliq(:,si));
            
        end
        
        % ASSIGN
        VARIABLES.SOIL.dwat = dwat;
        VARIABLES.SOIL.volliq = volliq;
        VARIABLES.SOIL.smp = smp;
        VARIABLES.SOIL.qlayer = qlayer;
        VARIABLES.SOIL.hor_drainage = hor_drainage;
        VARIABLES.SOIL.hor_drainage_lay = hor_drainage_lay;
        VARIABLES.SOIL.layeruptake = layeruptake;
        VARIABLES.SOIL.layeruptake_all = layeruptake_all;
        
        % RECOMPUTE MASS BALANCE INCLUDING THE FLUX BACK FROM
        % INFILTRATION SOLUTION
        
        for si=1:num_species
            if SWITCHES.litter(si)
                [VARIABLES] = FLUXES_WATER_SOIL_LITTER_BACK (VARIABLES, VERTSTRUC,...
                    PARAMS, CONSTANTS, si);
            else
                [VARIABLES] = FLUXES_WATER_SOIL_BACK (VARIABLES,...
                    VERTSTRUC, PARAMS, CONSTANTS, si);
            end
        end
        
        % assign
        qinfl = VARIABLES.SOIL.qinfl;
        qinflL = VARIABLES.SOIL.qinfl;
        net_qinflL = VARIABLES.SOIL.net_qinflL;
        drainlitter = VARIABLES.SOIL.drainlitter;
        volliqli = VARIABLES.SOIL.volliqli;
        
        if (SWITCHES.soilheat_on)
            
            
            % Soil Temperature Solution
            volice = 0;
            
            for si=1:num_species
                
                Ginto = VARIABLES.SOIL.Gs(si);
                
                [VARIABLES] = SOILHEAT (Ginto, VARIABLES, VERTSTRUC, PARAMS, CONSTANTS, si);
            end
            
            Ts = VARIABLES.SOIL.Ts;
            cpv = VARIABLES.SOIL.cpv;
        end
        
        % ************************************************************************
        %     [PARAMS, VARIABLES] = Nitrogen_Plant(PARAMS, FORCING, VARIABLES, CONSTANTS,nspecies);
        % ************************************************************************
        if SWITCHES.soilCN_on
            
            %             VARIABLES.Fertilized=1;
            
            %% RohitN Start
            if SWITCHES.NCanAlMod ==1
                %
                %                 td = 24 * 60 / CONSTANTS.timestep;
                %                 tq = tt + (CROP_GROWTH.DOY_start-1) * td;
                
                ph_type = PARAMS.Photosyn.ph_type;
                
                time_inter = 15;
                
                optimization=1;
                
                for si=1:num_species
                    
                    if (ph_type(si)~=1)
                        
                        VARIABLES.si=si;
                        
                        N_CAN_model;
                        
                        VARIABLES.CANOPY.Fert_req=Fert_req;
                        
                    end
                    
                end
                
            end
            
            
            % Fertilizer Selection
            
            Select_Fert;
            
            
            % For nitrogen remobilization
            if tt == 1
                STORAGE.UP_nit_m2_store = zeros(nl_soil,nl_soil,1);
                STORAGE.UP_amm_m2_store = zeros(nl_soil,nl_soil,1);
            end
            
            
            [VARIABLES, SWITCHES, PARAMS, FORCING] = ...
                core_N(rootfr, PARAMS, SWITCHES, VARIABLES, FORCING, CONSTANTS, VERTSTRUC, STORAGE);
            
            
            CN_STORE_DATA ();
            
        end
        
        if (SWITCHES.entropy_on)
            [SSresults] = ...
                COMPUENTROPY (SWcandir_in, SWcandir_out, SWcandif_in, SWcandif_out,...
                SWsoildir_in, SWsoildir_out, SWsoildif_in, SWsoildif_out,...
                SWout, fdiff,LWabs_canM, LWabs_soilM, LWemit_soil, LWemit_sun, LWemit_shade,...
                LWout, Tsurf, FORCING,SWITCHES, CONSTANTS, PARAMS,...
                VARIABLES, VERTSTRUC);
        end
        
        [VARIABLES] = MASS_BALANCE (VARIABLES, CONSTANTS, PARAMS, FORCING, SWITCHES, VERTSTRUC, tt);
        
        STORE_DATA;
        
    end
    
    NoDay_prev = NoDay_prev + NoDay;                     % Previous year NoDay
    
    toc
    
    % Timebar 3/3
    close(hh);
    
    SaveResFiles;
    
end


