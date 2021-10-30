
%     Initialize for Crop growth Model parameters
    CROP_GROWTH.LAI_ref=0.1*ones(24*NoDay,1);
    CROP_GROWTH.TotCmass_ref=0*ones(24*NoDay,1);
    CROP_GROWTH.Daily_SM_prev = 0;
    CROP_GROWTH.value = 0;
    CROP_GROWTH.lock=1;
    CROP_GROWTH.No_count=0;
    CROP_GROWTH.Vol_SM_all = null(24);
    CROP_GROWTH.doy1=zeros(1, num_species);
    CROP_GROWTH.doy2=zeros(1, num_species);
    CROP_GROWTH.doy3=zeros(1, num_species);
    CROP_GROWTH.doy4=zeros(1, num_species);
    CROP_GROWTH.doy5=zeros(1, num_species);
    CROP_GROWTH.ph_can=zeros(1, num_species);       %Photosynthesis
    VARIABLES.dim_temp=zeros(1, num_species);
    
    CROP_GROWTH.LeafMass_can = zeros(1, num_species);
    CROP_GROWTH.StemMass_can = zeros(1, num_species);
    CROP_GROWTH.RootMass_can = zeros(1, num_species);
    CROP_GROWTH.GrainMass_can = zeros(1, num_species);
    
    CROP_GROWTH.LeafMass_can_prev=zeros(1, num_species);
    CROP_GROWTH.RootMass_can_prev=zeros(1, num_species);
    CROP_GROWTH.StemMass_can_prev=zeros(1, num_species);
    CROP_GROWTH.GrainMass_can_prev=zeros(1, num_species);
    
    CROP_GROWTH.LAIsim_can = zeros(1, num_species);
    
    TotalMass_can = 0.0;
    AboveMass_can = 0.0;
    
    CROP_GROWTH.LeafBiomass = zeros(1, num_species);
    CROP_GROWTH.StemBiomass = zeros(1, num_species);
    CROP_GROWTH.RootBiomass = zeros(1, num_species);
    CROP_GROWTH.GrainBiomass = zeros(1, num_species);
    
    CROP_GROWTH.TotalBiomass_can = zeros(1, num_species);
    CROP_GROWTH.AboveBiomass_can = zeros(1, num_species);
    CROP_GROWTH.GrainYield = zeros(1, num_species);
    
     % Initialize parameters for Irrigation
    
    VARIABLES.CANOPY.TR_can = 0;
    VARIABLES.SOIL.Esoil = 0;
    VARIABLES.SOIL.ppt_ground = 0;
    VARIABLES.SOIL.ET_accu = 0;
    VARIABLES.CANOPY.psil_mean = 0;
    VARIABLES.CANOPY.gsv_mean = 0;
    VARIABLES.CANOPY.Tl_mean = Ta_in(1)*ones(1, num_species);
    VARIABLES.CANOPY.Tl_mean_up = Ta_in(1)*ones(1, num_species) + 2;
    VARIABLES.CANOPY.Tl_mean_lo = Ta_in(1)*ones(1, num_species);
    
    CROP_GROWTH.New_T = 0;
    CROP_GROWTH.time_delay=0;
    CROP_GROWTH.No_time=0;
    
%     Initialize the NCan variables
    N_Fert_DOY_lag = 0;
    Qabs_sun_NCan = [];
    Tl_sun_NCan =  [];
    Ci_sun_NCan =  [];
    Qabs_shade_NCan = [];
    Tl_shade_NCan = [];
    Ci_shade_NCan =  [];
    LAI_sun_can_NCan = [];
    LAI_shade_can_NCan = [];
    fLAIz_can_NCan= [];
    VARIABLES.NCan.K_ec=zeros(1, num_species);
    VARIABLES.NCan.K_rub_leaf=zeros(1, num_species);
    VARIABLES.NCan.K_lhc_leaf=zeros(1, num_species);
    VARIABLES.NCan.fn=zeros(1, num_species);
    
    VARIABLES.NCan.min_Vmax=8;
    VARIABLES.NCan.max_Vmax=65;
    VARIABLES.NCan.min_chl=0;
    VARIABLES.NCan.max_chl=1200;
    
    VARIABLES.NCan.N_can = zeros(1, num_species);
    VARIABLES.NCan.Vcmax_vert = PARAMS.Photosyn.Vmax_C4 .* ones (nl_can, num_species);
    VARIABLES.NCan.Vcmax_vert_x = PARAMS.Photosyn.Vmax_C4 .* ones (nl_can, num_species);
    VARIABLES.NCan.Chl = zeros(1, num_species);
    VARIABLES.NCan.Nleaf = zeros(nl_can, num_species);
    VARIABLES.NCan.Nrub = zeros (nl_can, num_species);
    VARIABLES.NCan.N_uptake_can = zeros(1, num_species);
    
    N_uptake_can =zeros(1, num_species);
    VARIABLES.NCan.Nall_start=zeros(1, num_species);
    CROP_GROWTH.N_Demand=zeros(1, num_species);
    VARIABLES.CANOPY.N_Demand = zeros(1, num_species);
    
    VARIABLES.NCan.Nl=zeros(1, num_species);
    VARIABLES.NCan.Nr=zeros(1, num_species);
    VARIABLES.NCan.Ns=zeros(1, num_species);
    VARIABLES.NCan.Ng=zeros(1, num_species);
    
    K_ec=0.8;
    K_rub_leaf=0.13;
    K_lhc_leaf=0.25;
    
    % IW=0;
    VARIABLES.SOIL.IW = zeros(1, num_species);
    
    save ('NCAN_opt_variables.mat', 'K_ec', 'K_rub_leaf', 'K_lhc_leaf');
    
    VARIABLES.NCan.K_ec=0.8*ones(1, num_species);
    VARIABLES.NCan.K_rub_leaf=0.13*ones(1, num_species);
    VARIABLES.NCan.K_lhc_leaf=0.25*ones(1, num_species);
    
%     Initialize Fertilizer parameters
    
    VARIABLES.UP_N_m2_sp=zeros(1, num_species);
    VARIABLES.CANOPY.N_req=0;
    VARIABLES.CANOPY.Net_Nreq=zeros(1, num_species);
    VARIABLES.CANOPY.Fert_req=0;
    
    Fert_req=0;
    
    VARIABLES.CANOPY.N_Fert_DOY=nan; % [DOY]
    VARIABLES.N_Fert_DOY_x=nan;
    VARIABLES.N_Fert_amm=nan; % [gN/m2]
    VARIABLES.N_Fert_nit=nan; % [gN/m2]
    VARIABLES.N_Fert_urea=nan; % [gN/m2]
    
    VARIABLES.CANOPY.DeathRate_prev = zeros(1, num_species);
    VARIABLES.CANOPY.RootDeathRate_prev = zeros(1, num_species);
    
    CROP_GROWTH.DeathRate = zeros(1, num_species);
    CROP_GROWTH.RootDeathRate = zeros(1, num_species);
    
    VARIABLES.CANOPY.DeathRate = zeros(1, num_species);
    VARIABLES.CANOPY.RootDeathRate = zeros(1, num_species);
    
    VARIABLES.Fert_Type=zeros(1, nspecies);
    
    %######################################%
    
    for si=1:num_species
        
        file_name_Py   = sprintf('YrwdSpec_%sY%s.mat', num2str(si), num2str(Run_years(yy)));
        
        load('.\Temps\temporary.mat', 'config_path');
        
        YrPath=fullfile(config_path, 'YearFr');
        
        fullpath_Py    =...
            fullfile(YrPath, file_name_Py);
        
        load(fullpath_Py, 'PlantingDOY', 'HarvestDOY', 'Vegetation')
        
        pinds = find(intdoy == PlantingDOY);
        hinds = find(intdoy == HarvestDOY);
        
        if SWITCHES.Onset==1
            CROP_GROWTH.plantingDate(si) = 0;
        else
            CROP_GROWTH.plantingDate(si) = intdoy(pinds(1))+ NoDay_prev;
        end
        
        CROP_GROWTH.harvestingDate(si) = intdoy(hinds(1))+ NoDay_prev;
        
        if strcmpi(Vegetation, 'Yes')
            SWITCHES.plants(si)=1;
            
        else
            SWITCHES.plants(si)=0;
            CROP_GROWTH.plantingDate(si)=0;
            CROP_GROWTH.harvestingDate(si)=0;
        end
        
        load(fullpath_Py, 'Irrigated', 'IrrType', 'IrrData')
        
        if Irrigated==0
            SWITCHES.Irr(si) = 0;
            VARIABLES.Irr.Type(si)=0;
        else
            SWITCHES.Irr(si) = 1;
            
            if strcmpi(IrrType, 'Manually')
                VARIABLES.Irr.Type(si)=1;
            elseif strcmpi(IrrType, 'LWP based')
                VARIABLES.Irr.Type(si)=2;
            elseif strcmpi(IrrType, 'CWSI based')
                VARIABLES.Irr.Type(si)=3;
            elseif strcmpi(IrrType, 'ET based')
                VARIABLES.Irr.Type(si)=4;
            end
        end
        
        if VARIABLES.Irr.Type(si)==1
            VARIABLES.Irr.DOY(si)=cell2mat(IrrData(:,1));
            VARIABLES.Irr.Amount(si)=cell2mat(IrrData(:,2));
        else
            VARIABLES.Irr.DOY(si)=0;
            VARIABLES.Irr.Amount(si)=0;
        end
        
        load(fullpath_Py, 'litter_depth')
        
        if litter_depth > 0
            SWITCHES.litter(si) = 1;          %  1 = Include litter dynamics. 0 does not include litter dynamics
        else
            SWITCHES.litter(si) = 0;          %  1 = Include litter dynamics. 0 does not include litter dynamics
        end
        
        if SWITCHES.litter(si)
            VARIABLES.SOIL.dzlit_m(si) = litter_depth;                            % Litter thickness in [m]
        else
            VARIABLES.SOIL.dzlit_m(si) = 0.;                                      % Litter thickness in [m]
        end
        
    end
    
    VARIABLES.SOIL.zliqsl = (VARIABLES.SOIL.dzlit_m*1000).*volliqliinit;
    VARIABLES.SOIL.zliqsl_prev = VARIABLES.SOIL.zliqsl;
    VARIABLES.SOIL.zsn = VARIABLES.SOIL.zliqsl;
    VARIABLES.SOIL.wliqsl = (VARIABLES.SOIL.zliqsl/1000).*PARAMS.Soil.rho_liq;
    VARIABLES.SOIL.wsn = VARIABLES.SOIL.wliqsl;
    