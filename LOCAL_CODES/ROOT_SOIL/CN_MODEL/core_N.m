% Dongkook Woo - Comment
% function []= core_N(VARIN)
% % Generate char for root cut type
% % Decode info for root cutting
% if VARIN(4) == 1
%    strcut = 'NN';
% elseif VARIN(4) == 2
%    strcut = 'OH';
% elseif VARIN(4) == 3
%    strcut = 'AH1';
% elseif VARIN(4) == 4
%    strcut = 'AH2';
% end
%
%
% % Generate VAR
% infile=[num2str(VARIN(1))  '_CUT' strcut '_HR' num2str(VARIN(2)) '_nsp' num2str(VARIN(3))];
% outfile1 = ['CN_' infile];
% outfile2 = ['CNSS_' infile];
%
% %clc;
% filename = ['./MLCAN_SA/' infile '.mat'];
% % Code Library Paths
%     addpath('./CN_MODEL/');
%
% if isempty(dir(filename))
%      disp(sprintf('file missing: %s\n',filename));
%      return;
% else
%     load (filename,'volliq_store','Ts_store','qlayer_store','wuptake_all_store',...
%     'wuptake_store','TR_can_all_store','smp_store','SWC1_in','SWC2_in','SWC3_in','Ts1_in',...
%     'Ts2_in','Ts3_in','Ts4_in','Ts5_in','Ts6_in','Ts7_in','Ts8_in','Ts9_in',...
%     'Ts10_in','dzsmm','zhsmm','znsmm','nl_soil','theta_dry',...
%     'porsl', 'rootfr','CONSTANTS','PARAMS','VERTSTRUC',...
%     'volliqli_store','Tli_store','qinflL_store','psili_store','Ts_store');
%     load ('tsoil_corr.mat','fTdl');
% end
% Dongkook Woo - Comment end

% Dongkook Woo insert
function  [VARIABLES, SWITCHES, PARAMS, FORCING]  = ...
    core_N(rootfr, PARAMS, SWITCHES, VARIABLES, FORCING, CONSTANTS, VERTSTRUC, STORAGE)
%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                              FUNCTION CODE                            %%
%%               Loading Parameters and Forcings for CN model            %%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Juan Quijano                                           %
%   Modified by  : Dongkook Woo                                           %
%   Date         : January 11, 2014                                       %
%-------------------------------------------------------------------------%
%                                                                         %
%   Loading parameters and forcings Soil Carbon and Nitrogen Model        %
%                                                                         %
%-------------------------------------------------------------------------%
%%
if VARIABLES.timestep == 1
    % Dongkook Woo insert end
    load './Temps/temp_model.mat'...
        'working_forcings' 'Soil_C_pool1'
    load (working_forcings)
    
    nl_soil=PARAMS.nl_soil;
    
    nspecies=PARAMS.CanStruc.nspecies;
    
    N=length(CNratio_species1_crop);
    
    %VARIABLES.SOIL.litterthickness = 0.03;
    VARIABLES.SOIL.litterthickness = VARIABLES.SOIL.dzlit_m;
    %N = 17520;
    
    %*************************************************************************%
    %                             SWITCHES                                   *%
    %*************************************************************************%
    % switches ();
    % SWITCHES.CN_type = 0;             % 1 = 3 pools Include humus
    % % 2 = 3 pools steady stated in humus
    % % 0 = 2 pools humus and litter are coupled
    if Soil_C_pool1 == 1
        SWITCHES.CN_type = 0;   % 1 = 3 pools Include humus
        % 2 = 3 pools steady stated in humus
        % 0 = 2 pools humus and litter are coupled
    elseif Soil_C_pool1 == 0
        SWITCHES.CN_type = 1;
    end
    
    %SWITCHES.initialcond = 3;         % 1 = load data with file 2. fill with desire numbers 3. same as Cl_bloi
    SWITCHES.CN.Bioturbation = 1;     % 0 No Bioturbation, 1 Including Bioturbation (include litter)
    %     SWITCHES.Recycling = 0;           % 1 = yes, there is recycling of NO3 in the soil by HR
    % 0 = No, there is not recycling of NO3 in the soil by HR
    %SWITCHES.CN.consmethod = 2;       % Method to calulate bioturbation in the computation of constants kl,kd
    % 1. using bioturbation fluxes
    % 2. using an exponential function
    %SWITCHES.comconstants = 0;        % 1. Compute Constants. 0. Does not compute Constants (load from file)
    
    % SWITCHES.save = 0;
    % SWITCHES.savetype = 0;
    % SWITCHES.fitlitter = 0;
    load './Temps/temp_model.mat'...
        'N_denit' 'N_Adepo' 'N_Uptake_RB' 'N_Fix' 'N_Remo'
    SWITCHES.CN.Denitrification = N_denit; % 1 = Including denitrification 0= No denitrification
    SWITCHES.CN.AtmosDeposition = N_Adepo; % 1 = Including Atmospheric Nitrogen Deposition 0= No Atmospheric Nitrogen Deposition
    
    SWITCHES.CN.NupRootBiomass  = N_Uptake_RB; % 1 = Computing N uptake with root biomass
    % 0 = Computing N uptake with plant C dynamics
    SWITCHES.CN.leach_on        = 1;%       leach_on = Leaching ON or OF [1 0]
    SWITCHES.CN.N_Fix           = N_Fix; % 1= including fixation 0 no
    SWITCHES.CN.N_Remo          = N_Remo; % 1=including fixation 0 no
    SWITCHES.Active_Nuptake     = 1;
    SWITCHES.Recycling          = 0; % Switch that decides if recycling of nitrate is on or off
    %*************************************************************************%
    %                            PARAMETERS                                  *%
    %*************************************************************************%
    
    % FOR CARBON NITROGEN SOIL MODEL
    %PARAMS.CN.factorrec = 0.3;     %
    %Fraction of C and N in leaves and fine roots
    %PARAMS.CN.fCleaves=0.55;
    %PARAMS.CN.fCroots=0.5;
    %PARAMS.CN.Nrecyc=0.0;      % Amount of Nitrogen that is recycled before leaves drop
    %PARAMS.CN.CR_ratio=0.36;   % Portion of C in roots compared to leaves
    
    load './Temps/temp_model.mat'...
        'DOY_start' 'DOY_end'
    PARAMS.CN.DOY_start = DOY_start;
    PARAMS.CN.DOY_end = DOY_end;
    
    %  MODIFICATIONS TO INCLUDE IN BIG CODE
    
    filename_P= ['ParSpec_', num2str(1), '.mat'];
        
    load('.\Temps\temporary.mat', 'config_path');
    
    PrPath=fullfile(config_path, 'Par');
    
    fullpath_P= fullfile(PrPath, filename_P);
    
    load(fullpath_P, 'para_soilCN');
    
    
    if isempty(cell2mat(para_soilCN(8,2)))
        PARAMS.CN.a_Amm = cell2mat(para_soilCN(8,3)); % [cm2/year]
    else
        PARAMS.CN.a_Amm = cell2mat(para_soilCN(8,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(9,2)))
        PARAMS.CN.a_Nit = cell2mat(para_soilCN(9,3)); % [cm2/year]
    else
        PARAMS.CN.a_Nit = cell2mat(para_soilCN(9,2)); % [cm2/year]
    end
    %PARAMS.CN.DEMp = 0.2;         % [gN / m^3 / d]
    %PARAMS.CN.DEMm = 0.5;         % [gN / m^3 / d]
    
    
    %PARAMS.CN.kn = 0.6/100;       % [m^3 / d / gC]
    
    %    PARAMS.CN.ki_Amm = 1;         % [m^3 / d / gC]
    %    PARAMS.CN.ki_Nit = 1;         % [m^3 / d / gC]
    %    PARAMS.CN.ku_Amm = 0.1;       % Maximum Immobilization [m^3 / d / gC]
    %    PARAMS.CN.ku_Nit = 0.1;       % Maximum Immobilization [m^3 / d / gC]
    %    PARAMS.CN.rhmin = 0.2;
    
    if isempty(cell2mat(para_soilCN(4,2)))
        PARAMS.CN.ki_Amm = cell2mat(para_soilCN(4,3)); % [cm2/year]
    else
        PARAMS.CN.ki_Amm = cell2mat(para_soilCN(4,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(5,2)))
        PARAMS.CN.ki_Nit = cell2mat(para_soilCN(5,3)); % [cm2/year]
    else
        PARAMS.CN.ki_Nit = cell2mat(para_soilCN(5,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(6,2)))
        PARAMS.CN.ku_Amm = cell2mat(para_soilCN(6,3)); % [cm2/year]
    else
        PARAMS.CN.ku_Amm = cell2mat(para_soilCN(6,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(7,2)))
        PARAMS.CN.ku_Nit = cell2mat(para_soilCN(7,3)); % [cm2/year]
    else
        PARAMS.CN.ku_Nit = cell2mat(para_soilCN(7,2)); % [cm2/year]
    end
    if isempty(cell2mat(para_soilCN(14,2)))
        PARAMS.CN.rhmin = cell2mat(para_soilCN(14,3)); % [cm2/year]
    else
        PARAMS.CN.rhmin = cell2mat(para_soilCN(14,2)); % [cm2/year]
    end
    
    %PARAMS.CN.rr = 0.7;
    PARAMS.CN.koae = 1 ; % Organic Assimilation Efficiency parameter
    %PARAMS.CN.Ts_reduc_on = 0;
    %PARAMS.CN.Ts_max = 10;
    %PARAMS.CN.leach_on = 1;
    %     if (SWITCHES.soilCN_on && SWITCHES.CN.Bioturbation)
    %         CNb(1:PARAMS.nl_soil+1) = 11.8;
    %         load ('tsoil_corr.mat','fTdl');
    %         PARAMS.CN.fTdl = fTdl;
    %     else
    %         CNb(1:PARAMS.nl_soil) = 11.8;
    %     end
    %     PARAMS.CN.CNb= CNb(:);
    %PARAMS.CN.CNh = 22;
    
    
    %bioturbation
    %PARAMS.CN.bio.biofactor = 0.5; % Factor used to compute Bioturbation rate[gr/m2]
    % Is a factor of the Total Biomass that comes
    % into the soil
    %PARAMS.kbio = 0.0000215;
    %PARAMS.CN.bio.dsoil = 2000*1000;%         % Bulk density Mineral Soil [gr/m3]
    %PARAMS.CN.bio.lini = 2;%         % layer where difussion starts
    %PARAMS.CN.bio.lfin = 12;%         % layer where difussion finish
    %PARAMS.CN.Clitter = 60000;        % [gr/m3]
    %PARAMS.CN.layerbio = 7;   % Maximum layer for bioturbation
    %PARAMS.CN.D = 4;         % [cm2/year]
    
    %PARAMS.CN.leafspan = 3;               % leaf span shrubs [years]
    
    %PARAMS.K2=0.02;
    
    if isempty(cell2mat(para_soilCN(12,2)))
        PARAMS.CN.K2 = cell2mat(para_soilCN(12,3)); % [cm2/year]
    else
        PARAMS.CN.K2 = cell2mat(para_soilCN(12,2)); % [cm2/year]
    end
    
    
    load './Temps/temp_model.mat'...
        'opt_root_para' 'N_Uptake_RB' 'Start_Y'
    
    
    if isempty(cell2mat(para_soilCN(10,2)))
        PARAMS.CN.D = cell2mat(para_soilCN(10,3)); % [cm2/year]
    else
        PARAMS.CN.D = cell2mat(para_soilCN(10,2)); % [cm2/year]
    end
    if N_Uptake_RB == 1
        PARAMS.Soil.rootarq.diam = cell2mat(opt_root_para(1,2:end)); % Root Diameter [mm]
        PARAMS.Soil.rootarq.dist = cell2mat(opt_root_para(2,2:end)); % Distribution
        PARAMS.Soil.rootarq.SRL  = cell2mat(opt_root_para(3,2:end));        % Specific Root Length [m/gr]
        if isempty(cell2mat(para_soilCN(19,2)))
            PARAMS.Soil.rootarq.deltarZD = cell2mat(para_soilCN(17,3)); % Definition zone of disturbance [mm]
        else
            PARAMS.Soil.rootarq.deltarZD = cell2mat(para_soilCN(17,3)); % Definition zone of disturbance [mm]
        end
    elseif N_Uptake_RB == 0
        PARAMS.Soil.rootarq.diam = nan;         % Root Diameter [mm]
        PARAMS.Soil.rootarq.dist = nan;         % Distribution
        PARAMS.Soil.rootarq.SRL = nan;          % Specific Root Length [m/gr]
        PARAMS.Soil.rootarq.deltarZD = nan;     % Definition zone of disturbance [mm]
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%--Needed to change--%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load ('.\LOCAL_CODES\ROOT_SOIL\CN_MODEL\jaeger.mat');
    %PARAMS.CN.jaeger = jaeger_1942;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % PARAMS.vmaxnit = 0.5;%0.037; %[mumolNO3/g_rr/hr]
    % PARAMS.kmnit = 72;%15; %[mmolNO3/m^3]
    % PARAMS.vmaxamm = 0.053; %[mumolNO3/g_rr/hr]
    % PARAMS.kmamm = 15; %[mmolNO3/m^3]
    
    if isempty(cell2mat(para_soilCN(18,2)))
        PARAMS.vmaxnit = cell2mat(para_soilCN(18,3)); %[mumolNO3/g_rr/hr]
    else
        PARAMS.vmaxnit = cell2mat(para_soilCN(18,2)); %[mumolNO3/g_rr/hr]
    end
    if isempty(cell2mat(para_soilCN(19,2)))
        PARAMS.kmnit = cell2mat(para_soilCN(19,3)); %[mmolNO3/m^3]
    else
        PARAMS.kmnit = cell2mat(para_soilCN(19,2)); %[mmolNO3/m^3]
    end
    if isempty(cell2mat(para_soilCN(20,2)))
        PARAMS.vmaxamm = cell2mat(para_soilCN(20,3)); %[mumolNO3/g_rr/hr]
    else
        PARAMS.vmaxamm = cell2mat(para_soilCN(20,2)); %[mumolNO3/g_rr/hr]
    end
    if isempty(cell2mat(para_soilCN(21,2)))
        PARAMS.kmamm = cell2mat(para_soilCN(21,3)); %[mmolNO3/m^3]
    else
        PARAMS.kmamm = cell2mat(para_soilCN(21,2)); %[mmolNO3/m^3]
    end
    
    if isempty(cell2mat(para_soilCN(2,2)))
        PARAMS.CN.kn = cell2mat(para_soilCN(2,3));
    else
        PARAMS.CN.kn = cell2mat(para_soilCN(2,2));
    end
    if isempty(cell2mat(para_soilCN(1,2)))
        PARAMS.CN.CNb = ones(PARAMS.nl_soil+1,1)*cell2mat(para_soilCN(1,3));
    else
        PARAMS.CN.CNb = ones(PARAMS.nl_soil+1,1)*cell2mat(para_soilCN(1,2));
    end
    
    % PARAMS.kbio = 0.0000215;
    % PARAMS.CN.layerbio = 7;   % Maximum layer for bioturbation
    if isempty(cell2mat(para_soilCN(11,2)))
        PARAMS.kbio = cell2mat(para_soilCN(11,3));
    else
        PARAMS.kbio = cell2mat(para_soilCN(11,2));
    end
    
    %CN_nutrients_parameters;
    %TBMla = FORCING.TBMla;
    %CNh = PARAMS.CN.CNh;
    if isempty(cell2mat(para_soilCN(13,2)))
        PARAMS.CN.CNh = cell2mat(para_soilCN(13,3));
    else
        PARAMS.CN.CNh = cell2mat(para_soilCN(15,2));
    end
    if isempty(cell2mat(para_soilCN(3,2)))
        PARAMS.CN.rr = cell2mat(para_soilCN(3,3));
    else
        PARAMS.CN.rr = cell2mat(para_soilCN(3,2));
    end
    
    
    nspecies = PARAMS.CanStruc.nspecies;
    
    
    for si=1:nspecies
        filename_P= ['ParSpec_', num2str(si), '.mat'];
        
        load('.\Temps\temporary.mat', 'config_path');
        
        PrPath=fullfile(config_path, 'Par');
        
        fullpath_P= fullfile(PrPath, filename_P);
        
        load (fullpath_P, 'para_soilCN_crop', 'para_soilCN_crop_s')
        
        file_name_Py   = sprintf('YrwdSpec_%sY%s.mat', num2str(si), num2str(Start_Y));
        
        load('.\Temps\temporary.mat', 'config_path');
        
        YrPath=fullfile(config_path, 'YearFr');
        
        fullpath_Py    =...
            fullfile(YrPath, file_name_Py);
        
%         load(fullpath_Py, 'num_root')
%         
        
        if N_Uptake_RB == 1
            
            if isempty(cell2mat(para_soilCN_crop(1,2)))
                PARAMS.CN.CR_ratio(si) = cell2mat(para_soilCN_crop(1,3));
            else
                PARAMS.CN.CR_ratio(si) = cell2mat(para_soilCN_crop(1,2));
            end
            
        elseif N_Uptake_RB == 0
            
            if isempty(cell2mat(para_soilCN_crop_s(1,2)))
                PARAMS.CN.fDEMamm(si) = cell2mat(para_soilCN_crop_s(1,3));
            else
                PARAMS.CN.fDEMamm(si) = cell2mat(para_soilCN_crop_s(1,2));
            end
            
            if isempty(cell2mat(para_soilCN_crop_s(2,2)))
                PARAMS.CN.fDEMnit(si) = cell2mat(para_soilCN_crop_s(2,3));
            else
                PARAMS.CN.fDEMnit(si) = cell2mat(para_soilCN_crop_s(2,2));
            end
            
            if isempty(cell2mat(para_soilCN_crop_s(3,2)))
                PARAMS.CN.fTortuo(si) = cell2mat(para_soilCN_crop_s(3,3));
            else
                PARAMS.CN.fTortuo(si) = cell2mat(para_soilCN_crop_s(3,2));
            end
            
            if isempty(cell2mat(para_soilCN_crop_s(4,2)))
                PARAMS.CN.fRd(si) = cell2mat(para_soilCN_crop_s(4,3));
            else
                PARAMS.CN.fRd(si) = cell2mat(para_soilCN_crop_s(4,2));
            end
            
            
            if isempty(cell2mat(para_soilCN_crop_s(5,2)))
                PARAMS.CN.NfI(si) = cell2mat(para_soilCN_crop_s(5,3));
            else
                PARAMS.CN.NfI(si) = cell2mat(para_soilCN_crop_s(5,2));
            end
            
            
            if isempty(cell2mat(para_soilCN_crop_s(6,2)))
                PARAMS.CN.NrI(si) = cell2mat(para_soilCN_crop_s(6,3));
            else
                PARAMS.CN.NrI(si) = cell2mat(para_soilCN_crop_s(6,2));
            end
            
        end
        
    end
    
    
    %fc= 0.39;     %Saxton et al (1986)
    %fHorO = 0.;
    %load constants.mat;
    PARAMS.CN.fHorO =0;
    load './Temps/temp_model.mat'...
        'dat_decom' 'dat_decom_litter'
    % kd=[kd(1); kd(:)];
    % PARAMS.CN.kd = kd(:);
    % kl=[kl(1); kl(:)];
    % PARAMS.CN.kl = kl(:);
    % kh=[kh(1); kh(:)];
    % PARAMS.CN.kh = kh(:);
    
    if SWITCHES.CN.Bioturbation==1
        
        kd=[dat_decom_litter(1,2); dat_decom(:,2)];
        PARAMS.CN.kd = kd(:);
        kl=[dat_decom_litter(1,3); dat_decom(:,3)];
        PARAMS.CN.kl = kl(:);
        kh=[dat_decom_litter(1,4); dat_decom(:,4)];
        PARAMS.CN.kh = kh(:);
    else
        kd=dat_decom(:,2);
        PARAMS.CN.kd = kd(:);
        kl=dat_decom(:,3);
        PARAMS.CN.kl = kl(:);
        kh=dat_decom(:,4);
        PARAMS.CN.kh = kh(:);
    end
    
    %     PARAMS.CN.kd = ones(12+1,1)*8.5*10^-3;     % [d^-1]
    %     PARAMS.CN.kl = ones(12+1,1)*6.5*10^-6;     % [m^3 / d / gC]
    %     PARAMS.CN.kh = ones(12+1,1)*2.5*10^-7;     % [m^3 / d / g?]
    
    
    load './Temps/temp_model.mat'...
        'N_Adepo' 'N_Adepo_amm' 'N_Adepo_nit'
    if N_Adepo == 1
        PARAMS.CN.N_Adepo_amm = N_Adepo_amm; % [g/m2/yr]
        PARAMS.CN.N_Adepo_nit = N_Adepo_nit; % [g/m2/yr]
    elseif N_Adepo == 0
        PARAMS.CN.N_Adepo_amm = nan;
        PARAMS.CN.N_Adepo_nit = nan;
    end
    
    %     Fertilizers
    
    if isempty(cell2mat(para_soilCN(15,2)))
        PARAMS.CN.VolUrea = cell2mat(para_soilCN(15,3));
    else
        PARAMS.CN.VolUrea = cell2mat(para_soilCN(15,2));
    end
    if isempty(cell2mat(para_soilCN(16,2)))
        PARAMS.CN.VolAmm = cell2mat(para_soilCN(16,3));
    else
        PARAMS.CN.VolAmm = cell2mat(para_soilCN(16,2));
    end
    
    %*************************************************************************%
    %                             Forcings                                   *%
    %*************************************************************************%
    
    %********************************************************************
    %  JUST FOR THIS CODE
    load './Temps/temp_model.mat'...
        'working_forcings'
    load (working_forcings)
    
    % Organize the data for simulation
    load './Temps/temp_model.mat'...
        'working_forcings' 'DOY_start' 'DOY_end'
    load(working_forcings,...
        'year_crop','doy_crop');
    
    years = unique(year_crop)';
    
    load('./Temps/temp_model.mat', 'Start_Y', 'End_Y')

    Run_years=Start_Y:End_Y;
    
    doys = [DOY_start:DOY_end];
    inds = find(ismember(year_crop, Run_years) & ismember(floor(round((doy_crop+1).*10^10)./10^10), doys));
    FORCING.doy_all=doy_crop(inds);
    
    if N_Uptake_RB == 1
        
        PARAMS.CN.CNveg{1} = CNratio_species1_crop(inds);
        PARAMS.CN.CNveg{2} = CNratio_species2_crop(inds);
        PARAMS.CN.CNveg{3} = CNratio_species3_crop(inds);
        PARAMS.CN.CNveg{4} = CNratio_species4_crop(inds);
        FORCING.TBMla(1,:) = BMla_species1_crop(inds)';
        FORCING.TBMla(2,:) = BMla_species2_crop(inds)';
        FORCING.TBMla(3,:) = BMla_species2_crop(inds)';
        FORCING.TBMla(4,:) = BMla_species2_crop(inds)';
        FORCING.LMI(1,:)   = LMI_species1_crop(inds);
        FORCING.LMI(2,:)   = LMI_species2_crop(inds);
        FORCING.LMI(3,:)   = LMI_species3_crop(inds);
        FORCING.LMI(4,:)   = LMI_species4_crop(inds);
        
        
        
    elseif N_Uptake_RB == 0
        
        FORCING.BMla(1,:)    = BMlaDEM_species1_crop(inds);
        FORCING.BMrb(1,:)    = BMrbDEM_species1_crop(inds);
        FORCING.CNabove(1,:) = CNratioA_species1_crop(inds);
        FORCING.CNbelow(1,:) = CNratioB_species1_crop(inds);
        FORCING.CNdem(1,:)   = NDEM_species1_crop(inds);
        
        FORCING.BMla(2,:)    = BMlaDEM_species2_crop(inds);
        FORCING.BMrb(2,:)    = BMrbDEM_species2_crop(inds);
        FORCING.CNabove(2,:) = CNratioA_species2_crop(inds);
        FORCING.CNbelow(2,:) = CNratioB_species2_crop(inds);
        FORCING.CNdem(2,:)   = NDEM_species2_crop(inds);
        
        FORCING.BMla(3,:)    = BMlaDEM_species3_crop(inds);
        FORCING.BMrb(3,:)    = BMrbDEM_species3_crop(inds);
        FORCING.CNabove(3,:) = CNratioA_species3_crop(inds);
        FORCING.CNbelow(3,:) = CNratioB_species3_crop(inds);
        FORCING.CNdem(3,:)   = NDEM_species3_crop(inds);
        
        FORCING.BMla(4,:)    = BMlaDEM_species4_crop(inds);
        FORCING.BMrb(4,:)    = BMrbDEM_species4_crop(inds);
        FORCING.CNabove(4,:) = CNratioA_species4_crop(inds);
        FORCING.CNbelow(4,:) = CNratioB_species4_crop(inds);
        FORCING.CNdem(4,:)   = NDEM_species4_crop(inds);
        
        
    end
    
    if VARIABLES.yy==1
        
        %*************************************************************************%
        %                        Initial Condition                               *%
        %*************************************************************************%
        
        % INITIALIZE NUTRIENTS. POOLS SIZE AND REQUIRED VARIABLES.
        VARIABLES.dCl_dT = 0;
        VARIABLES.DECl = 0;
        VARIABLES.BD = 0;
        
        load './Temps/temp_model.mat'...
            'Soil_C_pool1' 'nutrient_int' 'nutrient_int_litter'
        
        if SWITCHES.CN.Bioturbation==1
            
            if Soil_C_pool1 == 1
                Cl =[nutrient_int_litter(2) ; nutrient_int(:,2)];
                VARIABLES.Cl = Cl(:).*ones(1,nspecies);
                
                PARAMS.CN.Clitter=nutrient_int_litter(2).*ones(1,nspecies);
                
                Cb =[nutrient_int_litter(3) ; nutrient_int(:,3)];
                VARIABLES.Cb = Cb(:).*ones(1,nspecies);
                
                Ch =ones(PARAMS.nl_soil+1,1)*nan;
                VARIABLES.Ch =Ch(:).*ones(1,nspecies);
                
                VARIABLES.Amm = [nutrient_int_litter(4) ; nutrient_int(:,4)].*ones(1,nspecies);
                VARIABLES.Ammdd = [nutrient_int_litter(4) ; nutrient_int(:,4)].*ones(1,nspecies);
                
                VARIABLES.Nit = [nutrient_int_litter(5) ; nutrient_int(:,5)].*ones(1,nspecies);
                VARIABLES.Nitdd = [nutrient_int_litter(5) ; nutrient_int(:,5)].*ones(1,nspecies);
                
                VARIABLES.CNl = [nutrient_int_litter(6) ; nutrient_int(:,6)].*ones(1,nspecies);  % Initialize with CN min
                VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;
                
            elseif Soil_C_pool1 ==0
                Cl =[nutrient_int_litter(2) ; nutrient_int(:,2)];
                VARIABLES.Cl = Cl(:).*ones(1,nspecies);
                
                PARAMS.CN.Clitter=nutrient_int_litter(2).*ones(1,nspecies);
                
                Cb =[nutrient_int_litter(4) ; nutrient_int(:,4)];
                VARIABLES.Cb = Cb(:).*ones(1,nspecies);
                
                Ch =[nutrient_int_litter(3) ; nutrient_int(:,3)];
                VARIABLES.Ch =Ch(:).*ones(1,nspecies);
                
                VARIABLES.Amm = [nutrient_int_litter(5) ; nutrient_int(:,5)].*ones(1,nspecies);
                VARIABLES.Ammdd = [nutrient_int_litter(5) ; nutrient_int(:,5)].*ones(1,nspecies);
                
                VARIABLES.Nit = [nutrient_int_litter(6) ; nutrient_int(:,6)].*ones(1,nspecies);
                VARIABLES.Nitdd = [nutrient_int_litter(6) ; nutrient_int(:,6)].*ones(1,nspecies);
                
                VARIABLES.CNl = [nutrient_int_litter(7) ; nutrient_int(:,7)].*ones(1,nspecies);  % Initialize with CN min
                VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;
                
            end
            
        else
            if Soil_C_pool1 == 1
                Cl =nutrient_int(:,2);
                VARIABLES.Cl = Cl(:).*ones(1,nspecies);
                
                Cb =nutrient_int(:,3);
                VARIABLES.Cb = Cb(:).*ones(1,nspecies);
                
                Ch =ones(PARAMS.nl_soil,1)*nan;
                VARIABLES.Ch =Ch(:).*ones(1,nspecies);
                
                VARIABLES.Amm = nutrient_int(:,4).*ones(1,nspecies);
                VARIABLES.Ammdd = nutrient_int(:,4).*ones(1,nspecies);
                
                VARIABLES.Nit = nutrient_int(:,5).*ones(1,nspecies);
                VARIABLES.Nitdd = nutrient_int(:,5).*ones(1,nspecies);
                
                VARIABLES.CNl = nutrient_int(:,6).*ones(1,nspecies);  % Initialize with CN min
                VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;
                
            elseif Soil_C_pool1 ==0
                Cl = nutrient_int(:,2);
                VARIABLES.Cl = Cl(:).*ones(1,nspecies);
                
                Cb =nutrient_int(:,4);
                VARIABLES.Cb = Cb(:).*ones(1,nspecies);
                
                Ch = nutrient_int(:,3);
                VARIABLES.Ch =Ch(:).*ones(1,nspecies);
                
                VARIABLES.Amm =  nutrient_int(:,5).*ones(1,nspecies);
                VARIABLES.Ammdd = nutrient_int(:,5).*ones(1,nspecies);
                
                VARIABLES.Nit = nutrient_int(:,6).*ones(1,nspecies);
                VARIABLES.Nitdd = nutrient_int(:,6).*ones(1,nspecies);
                
                VARIABLES.CNl = nutrient_int(:,7).*ones(1,nspecies);  % Initialize with CN min
                VARIABLES.Nl = VARIABLES.Cl./VARIABLES.CNl;
                
            end
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%--Needed to change--%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try to fix it! -> ok
    %VARIABLES.SOIL.psil=VARIABLES.SOIL.smp(1);
    %VARIABLES.SOIL.psil;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %*************************************************************************%
    %                               Allocaition                              *%
    %*************************************************************************%
    addpath('./LOCAL_CODES/ROOT_SOIL/CN_MODEL');
    ALLOCATE ();
    
end


% load './Temps/temp_model.mat'...
%         'N_Fert' 'N_Fert_DOY' 'N_Fert_amm' 'N_Fert_nit' 'N_Fert_urea'
%

SWITCHES.CN.Fertilizer = VARIABLES.Fertilized; % 1 = Including Fertilizer application 0= No Fertilizer application              = 0;

nspecies = PARAMS.CanStruc.nspecies;

for si=1:nspecies
    
    N_Fert=VARIABLES.Fertilized(si);
    
    Fert_Type=VARIABLES.Fert_Type(si);
    
    if N_Fert == 1 && (Fert_Type == 1 || Fert_Type == 2 || Fert_Type == 3)
        
        PARAMS.CN.N_Fert_DOY(si)  = VARIABLES.N_Fert_DOY_x(si); % [DOY]
        PARAMS.CN.N_Fert_amm(si)  = VARIABLES.N_Fert_amm(si); % [gN/m2]
        PARAMS.CN.N_Fert_nit(si)  = VARIABLES.N_Fert_nit(si); % [gN/m2]
        PARAMS.CN.N_Fert_urea(si) = VARIABLES.N_Fert_urea(si); % [gN/m2]
        
    else
        
        PARAMS.CN.N_Fert_DOY(si)  = nan; % [DOY]
        PARAMS.CN.N_Fert_amm(si)  = nan; % [gN/m2]
        PARAMS.CN.N_Fert_nit(si)  = nan; % [gN/m2]
        PARAMS.CN.N_Fert_urea(si) = nan; % [gN/m2]
        
    end
end



%*************************************************************************%
%             Calculate Dynamics of Soil Carbon and Nitrogen             *%
%*************************************************************************%
[VARIABLES] = ...
    CN_Cycle_Porporato2(rootfr, PARAMS, SWITCHES, VARIABLES, FORCING, CONSTANTS, VERTSTRUC, STORAGE);



