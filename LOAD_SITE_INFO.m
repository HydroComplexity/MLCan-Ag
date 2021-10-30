%************************** Number of Species ****************************%
load './Temps/temp_model.mat'...
    'num_species' 'Sim_species' 'Sim_species_con'
% num_species : Number of species
% Sim_species = 1: Multi species
% Sim_species_con: Which species you would like to run

% if Sim_species == 1
%     kspecies = num_species;
% else % for the case you would like to simulate a selected species
%     kspecies = 5;
% end

kspecies = num_species;

%***************************** Set Years *********************************%
load './Temps/temp_model.mat'...
    'working_forcings'
load(working_forcings,...
    'year_crop');
Each_year = unique(year_crop)';

%************************ Set Initial Condition **************************%
load './Temps/temp_model.mat'...
    'root_init' 'root_init_litter'
Tsinit=root_init(:,3);                          % Initial soil temperature ['C]
Tslint=root_init_litter(1,3);                   % Initial litter temperature ['C]
volliqinit=root_init(:,2);                      % Initial soil moisture [-]
volliqliinit=root_init_litter(1,2);             % Initial litter moisture [-]

%************ Site Independent Constant and Fixed Values *****************%
% CONVERSION FACTORS
CONSTANTS.umoltoWm2 = (2.17*10^5) / 10^6;       % Radiation conversion
CONSTANTS.Wm2toumol = 1/CONSTANTS.umoltoWm2;    % Radiation conversion
CONSTANTS.mmH2OtoMPa = 9.8066e-06;              % Pressure conversion

% PHYSICAL CONSTANTS
CONSTANTS.R = 8.314;                            % [J mol^-1 K^-1]
CONSTANTS.R_kJ = CONSTANTS.R/1000;              % [kJ mol^-1 K^-1]
CONSTANTS.Lv = 44000;                           % latent heat of vaporization [J/mol]
CONSTANTS.Lv_kg = 2260*1000;                    % Lv in [J / kg]
CONSTANTS.Lf = 6000;                            % latent heat of vaporization [J/mol]
CONSTANTS.Lf_kg = 3.337e5;                      % Lf in [J / kg]

CONSTANTS.cp_mol = 29.3;                        % specific heat of air at constant pressure [J/mol/K]
CONSTANTS.cp_JkgK = 1012;                       % [J/kg/K]
CONSTANTS.boltz = 5.6697 * 10^-8;               % Stefan-Boltzmann constant [W/m^2/K^4]
CONSTANTS.vonk = 0.41;                          % von Karman constant
% The K?m? constant is often used in turbulence modeling,
% for instance in boundary-layer meteorology to calculate fluxes
% of momentum, heat and moisture from the atmosphere to the land surface.
% It is considered to be a universal (? = 0.41).

CONSTANTS.rho_dry_air = 1.2923;                 % [kg / m^3]
CONSTANTS.grav = 9.8;                           % [m / s^2]

CONSTANTS.Vw = 18;                              % Mols to kg conversion

CONSTANTS.cpl_JkgK = 1800;                      % specific heat capacity of leaves[J/kg/K]
CONSTANTS.rho_leaf = 600;                       % Density of Leaves [kg / m^3]

CONSTANTS.fPAR = 0.5;                           % [-] fraction incoming shortwave as PAR

% ROOT CUT INFORMATION
% wilting point
CONSTANTS.wilpoint = -1.5;                      % Wilting Point 1.5 MPa
% embolism
CONSTANTS.PLCa = 2.146;                         % Parameter a embolism equation
CONSTANTS.PLCb = -1.238;                        % Parameter b embolism equation



load './Temps/temp_model.mat'...
    'days_step' 'mins_step' 'hours_step'
CONSTANTS.timestep = days_step*24*60+hours_step*60+mins_step; % [minutes]
CONSTANTS.dtime = CONSTANTS.timestep*60;        % [s]




%% PARAMETERS

%*************************************************************************%
%                               Canopy                                    %
%*************************************************************************%
%********************* Independent of vegetation *************************%
% CANOPY STRUCTURE

filename_P= ['ParSpec_', num2str(1), '.mat'];

load('.\Temps\temporary.mat', 'config_path');

PrPath=fullfile(config_path, 'Par');

fullpath_Pi= fullfile(PrPath, filename_P);

load (fullpath_Pi,...
    'para_canopy_crop_fixed')
if isempty(cell2mat(para_canopy_crop_fixed(3,2)))
    PARAMS.CanStruc.hcan = cell2mat(para_canopy_crop_fixed(3,3));     % canopy height [m]
else
    PARAMS.CanStruc.hcan = cell2mat(para_canopy_crop_fixed(3,2));
end
if isempty(cell2mat(para_canopy_crop_fixed(4,2)))
    PARAMS.CanStruc.hobs = cell2mat(para_canopy_crop_fixed(4,3));     % observation height [m]
else
    PARAMS.CanStruc.hobs = cell2mat(para_canopy_crop_fixed(4,2));
end
if isempty(cell2mat(para_canopy_crop_fixed(5,2)))
    PARAMS.CanStruc.z0 = cell2mat(para_canopy_crop_fixed(5,3));       % canopy roughness length [m]
else
    PARAMS.CanStruc.z0 = cell2mat(para_canopy_crop_fixed(5,2));
end
PARAMS.CanStruc.d0 = 2/3 * PARAMS.CanStruc.hcan;                      % canopy displacement height [m]
PARAMS.CanStruc.VAratio = 3/1000;                                     % Ratio of volume to area [m]

% RADIATION
load (fullpath_Pi,...
    'para_radiation')
if isempty(cell2mat(para_radiation(1,2)))
    PARAMS.Rad.transmiss = cell2mat(para_radiation(1,3));             % atmospheric transmissivity
else
    PARAMS.Rad.transmiss = cell2mat(para_radiation(1,2));
end
if isempty(cell2mat(para_radiation(2,2)))
    PARAMS.Rad.epsv = cell2mat(para_radiation(2,3));                  % vegetation emissivity
else
    PARAMS.Rad.epsv = cell2mat(para_radiation(2,2));
end
if isempty(cell2mat(para_radiation(3,2)))
    PARAMS.Rad.epss = cell2mat(para_radiation(3,3));                  % soil emissivity
else
    PARAMS.Rad.epss = cell2mat(para_radiation(3,2));
end
if isempty(cell2mat(para_radiation(4,2)))
    PARAMS.Rad.epsa = cell2mat(para_radiation(4,3));                  % atmospheric emissivity
else
    PARAMS.Rad.epsa = cell2mat(para_radiation(4,2));
end
if isempty(cell2mat(para_radiation(5,2)))
    PARAMS.Rad.xx = cell2mat(para_radiation(5,3));                    % leaf angle dist param (Spherical is a good assumption (Shade and Golstein 2002)
else
    PARAMS.Rad.xx = cell2mat(para_radiation(5,2));
end
if isempty(cell2mat(para_radiation(6,2)))
    PARAMS.Rad.clump = cell2mat(para_radiation(6,3));                 % leaf clumping parameter
else
    PARAMS.Rad.clump = cell2mat(para_radiation(6,2));
end
if isempty(cell2mat(para_radiation(7,2)))
    PARAMS.Rad.Kdf = cell2mat(para_radiation(7,3));                   % extinction coeff for diffuse around 0.7 according to figure 15.2 (Campbell and Norman 1998)
else
    PARAMS.Rad.Kdf = cell2mat(para_radiation(7,2));
end

if isempty(cell2mat(para_radiation(8,2)))
    PARAMS.Rad.absorp_PAR = cell2mat(para_radiation(8,3));            % leaf absorptivity to PAR
else
    PARAMS.Rad.absorp_PAR = cell2mat(para_radiation(8,2));
end
if isempty(cell2mat(para_radiation(9,2)))
    PARAMS.Rad.absorp_NIR = cell2mat(para_radiation(9,3));            % leaf absorptivity to NIR
else
    PARAMS.Rad.absorp_NIR = cell2mat(para_radiation(9,2));
end
if isempty(cell2mat(para_radiation(10,2)))
    PARAMS.Rad.refl_PAR = cell2mat(para_radiation(10,3));             % PAR reflection coeff
else
    PARAMS.Rad.refl_PAR = cell2mat(para_radiation(10,2));
end
if isempty(cell2mat(para_radiation(11,2)))
    PARAMS.Rad.refl_NIR = cell2mat(para_radiation(11,3));             % NIR reflection coeff
else
    PARAMS.Rad.refl_NIR = cell2mat(para_radiation(11,2));
end
if isempty(cell2mat(para_radiation(12,2)))
    PARAMS.Rad.refl_soil = cell2mat(para_radiation(12,3));            % soil reflection coeff
else
    PARAMS.Rad.refl_soil = cell2mat(para_radiation(12,2));
end
PARAMS.Rad.trans_PAR = 1 - PARAMS.Rad.absorp_PAR - PARAMS.Rad.refl_PAR;
PARAMS.Rad.trans_NIR = 1 - PARAMS.Rad.absorp_NIR - PARAMS.Rad.refl_NIR;

if isempty(cell2mat(para_radiation(13,2)))
    PARAMS.Rad.refl_snow = cell2mat(para_radiation(13,3));            % Reflectivity of new snow.
else
    PARAMS.Rad.refl_snow = cell2mat(para_radiation(13,2));
end
if isempty(cell2mat(para_radiation(14,2)))
    PARAMS.Rad.refl_snow_old = cell2mat(para_radiation(14,3));        % How much reflectivity of snow decreases with age (per 12 days)
else
    PARAMS.Rad.refl_snow_old = cell2mat(para_radiation(14,2));
end

PARAMS.LWmethod = 0;                                                  % Method to compute in LW_ATTENUATION FUNCITON LW 0. Initial Darren. 1. Corrected
PARAMS.retainLW = 0;                                                  % Factor to compute LW, what proportion of LW radiation emmited by a canopy

% layer, remains there instead of leaving. Only if LWmethod =1
PARAMS.LWcom = 3;                                                     % Computation of longwave incoming 1.DATA 2.Boltzman 3.Boltzman Corrected


%% RESPIRATION
load (fullpath_Pi,...
    'para_respiration', 'para_microenvironment')
if isempty(cell2mat(para_respiration(1,2)))
    PARAMS.Resp.Ro = cell2mat(para_respiration(1,3));                 % [umol / m^2 / s]
else
    PARAMS.Resp.Ro = cell2mat(para_respiration(1,2));
end
if isempty(cell2mat(para_respiration(2,2)))
    PARAMS.Resp.Q10 = cell2mat(para_respiration(2,3));
else
    PARAMS.Resp.Q10 = cell2mat(para_respiration(2,2));
end

% Turbulence
if isempty(cell2mat(para_microenvironment(1,2)))
    PARAMS.MicroEnv.Cd = cell2mat(para_microenvironment(1,3));        % Drag coefficient [-]
else
    PARAMS.MicroEnv.Cd = cell2mat(para_microenvironment(1,2));
end

PARAMS.MicroEnv.alph = CONSTANTS.vonk/3;                              % mixing length parameter [-] / Katul et al (BLM, 2004, p. 84)


load (fullpath_Pi,...
    'para_photosynthesisC3_crop')
if isempty(cell2mat(para_photosynthesisC3_crop(5,2)))
    PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC3_crop(5,3)); % Vertical Distribution of Photosynthetic Capacity
else
    PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC3_crop(5,2));
end


%%
%*************************************************************************%
%                               SOIL                                      %
%*************************************************************************%

load('./Temps/temp_model', 'Start_Y');

file_name_Py   = sprintf('YrwdSpec_%sY%s.mat', num2str(1), num2str(Start_Y));

load('.\Temps\temporary.mat', 'config_path');

YrPath=fullfile(config_path, 'YearFr');

fullpath_Py    =...
    fullfile(YrPath, file_name_Py);

load (fullpath_Py,...
    'num_root', 'set_root_func')

if isstr(num_root) == 1
    PARAMS.nl_soil=str2num(num_root);
else
    PARAMS.nl_soil=num_root;
end

PARAMS.Soil.scaleza = cell2mat(set_root_func(1,2));
PARAMS.Soil.scalezb = cell2mat(set_root_func(2,2));

% SOIL PROPERTIES

load (fullpath_Pi,...
    'para_soil')

if isempty(cell2mat(para_soil(1,2)))
    VERTSTRUC.sand = cell2mat(para_soil(1,3));
else
    VERTSTRUC.sand = cell2mat(para_soil(1,2));
end
if isempty(cell2mat(para_soil(2,2)))
    VERTSTRUC.clay = cell2mat(para_soil(2,3));
else
    VERTSTRUC.clay = cell2mat(para_soil(2,2));
end
if isempty(cell2mat(para_soil(3,2)))
    PARAMS.Soil.Cd_soil = cell2mat(para_soil(3,3));                   % soil drag coefficient [-]
else
    PARAMS.Soil.Cd_soil = cell2mat(para_soil(3,2));
end
if isempty(cell2mat(para_soil(4,2)))
    PARAMS.Soil.z0 = cell2mat(para_soil(4,3));                        % soil surface roughness length [m]
else
    PARAMS.Soil.z0 = cell2mat(para_soil(4,2));
end

PARAMS.Soil.alphCN=0.5;                                               % Cranck Nicholson for Heat Solution


%% CONSTANTS
PARAMS.Soil.smpmin = -1.e8;                                           % restriction for min of soil poten. [mm]
PARAMS.Soil.wimp   = 0.05;                                            % water impremeable if porosity less than wimp [-]

PARAMS.Soil.scalek  = 0.5;                                            % length scale for the exponential decrease in Ksat [m]


% HEAT CAPACITIES [J / kg / K)]
% Dry Air
PARAMS.Soil.HC_air = 1.00464 * 10^3;
% Water
PARAMS.Soil.HC_liq = 4.188 * 10^3;
% Ice
PARAMS.Soil.HC_ice = 2.11727 * 10^3;
% DENSITIES [kg / m^3]
PARAMS.Soil.rho_liq = 1000;
PARAMS.Soil.rho_ice = 917;
% THERMAL CONDUCTIVITIES  [W / m / K]
PARAMS.Soil.TK_liq = 0.6;
PARAMS.Soil.TK_ice = 2.29;
PARAMS.Soil.TK_air = 0.023;
% FREEZING TEMP OF FRESH WATER [K]
PARAMS.Soil.Tf = 273.16;

% FOR LITTER SNOW
PARAMS.Tcr_atm = 273;                                                 % Atmospheric temperature above which snow
PARAMS.Soil.HC_snow = 0.1;                                            % Snow holding capacity of liquid water
PARAMS.Soil.rhosn_min = 50;                                           % Snow compaction parameter [kg/m3]
PARAMS.Soil.fcomp = 0.6;                                              %
PARAMS.Soil.rhod = 150.00;                                            % Snow compaction parameter [kg/m3]

PARAMS.Soil.c1 = -1.4e4;
PARAMS.Soil.thetals = 0.9;                                            % value of soil moisture litter at saturation
PARAMS.Soil.thetafc = 0.025;                                          % Value of litter soil moisture at field capacity
PARAMS.Soil.thetatr = 0.18;                                           % Value of soil moisture after which
% rl becomes negligible
PARAMS.Soil.psill = 35.3;                                             % parameter to compute psi litter
PARAMS.Soil.bl = 2.42;                                                % parameter to compute psi litter
PARAMS.Soil.bdl = 42.5;                                               % Bulk density of litter [kg/m3]
PARAMS.Soil.rhowater = 1000;                                          % Density of liquid water [1000 kg/m3]
PARAMS.Soil.TK_litter = 0.15;                                         % Litter Thermal Conductivity [W/m/k] or [J/s/m/k]
PARAMS.Soil.TD_litter = 5.7*10^(-7);                                  % Litter thermal diffusivity  [m/s]
PARAMS.Soil.slfactor = 0.1;                                           % Percentage above which is considered as snow for energy balance
PARAMS.Soil.VHC_litter = 0.3*10^(6);                                  % Volumetric Heat Capacity of Litter [J/m3/K]
PARAMS.Soil.km = 0.0001;%0.00004;                                     % parameter to compute drainage from litter
PARAMS.Soil.bm = 0.1;%2.3;                                            % parameter to compute drainage from litter
PARAMS.Soil.thetamin = 0.0001;                                        % value of soil moisture at which evaporation becomes negligible
PARAMS.Soil.ldif = 3.1*10^(-8);
PARAMS.Soil.sdif = 2.12*10^(-5);
PARAMS.Soil.kklitter = 1000;                                          % [1/cm] Radiation attenuation
PARAMS.Soil.kksnow = 0.8;                                             % [1/cm] Radiation attenuation

% ENTROPY COMPUTATION.
PARAMS.Entropy.c2 = 2.336;
PARAMS.Entropy.c3 = 0.26;

%*************************************************************************%
%********************* DEPENDENT OF VEGETATION TYPE **********************%
%*************************************************************************%

PARAMS.CanStruc.nspecies=num_species;


%**************************** Species 1 **********************************%

for ii=1:num_species
    
    % Leaf Parameters
    
    filename_P= ['ParSpec_', num2str(ii), '.mat'];
    
    load('.\Temps\temporary.mat', 'config_path');
    
    PrPath=fullfile(config_path, 'Par');
    
    fullpath_P= fullfile(PrPath, filename_P);
    
    load (fullpath_P,...
        'para_leaf_crop')
    
    if isempty(cell2mat(para_leaf_crop(1,2)))
        PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop(1,3));  % multiplicative factor for Fc and LE
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
    else
        if isstr(cell2mat(para_leaf_crop(1,2))) == 1
            PARAMS.CanStruc.LEfact(ii) = str2num(cell2mat(para_leaf_crop(1,2)));
        else
            PARAMS.CanStruc.LEfact(ii) = cell2mat(para_leaf_crop(1,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop(2,2)))
        PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop(2,3));   % multiplicative factor for H and LW
        % calculations (1 = fluxes from only
        % one side of leaf, 2 = fluxes from
        % both sides)
        
    else
        if isstr(cell2mat(para_leaf_crop(2,2))) == 1
            PARAMS.CanStruc.Hfact(ii) = str2num(cell2mat(para_leaf_crop(2,2)));
        else
            PARAMS.CanStruc.Hfact(ii) = cell2mat(para_leaf_crop(2,2));
        end
    end
    
    if isempty(cell2mat(para_leaf_crop(3,2)))
        PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop(3,3));
    else
        if isstr(cell2mat(para_leaf_crop(3,2))) == 1
            PARAMS.CanStruc.LWfact(ii) = str2num(cell2mat(para_leaf_crop(3,2)));
        else
            PARAMS.CanStruc.LWfact(ii) = cell2mat(para_leaf_crop(3,2));
        end
    end
    
    % CANOPY STRUCTURE
    if isempty(cell2mat(para_leaf_crop(4,2)))
        PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop(4,3)); % 1 = broad leaves, 2 = needles
    else
        if isstr(cell2mat(para_leaf_crop(4,2))) == 1
            PARAMS.CanStruc.leaftype(ii) = str2num(cell2mat(para_leaf_crop(4,2)));
        else
            PARAMS.CanStruc.leaftype(ii) = cell2mat(para_leaf_crop(4,2));
        end
    end
    
    load (fullpath_P,...
        'para_canopy_crop', 'para_leaf_crop_s')
    if isempty(cell2mat(para_canopy_crop(1,2)))
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop(1,3));     % leaf width or needle diameter [m]
    else
        PARAMS.CanStruc.ld(ii) = cell2mat(para_canopy_crop(1,2));
    end
    if isempty(cell2mat(para_canopy_crop(2,2)))
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop(2,3));     % shoot diameter for conifers
        % leaf width for broadleaved vegetation (= ld)
    else
        PARAMS.CanStruc.lw(ii) = cell2mat(para_canopy_crop(2,2));
    end
    
    % PARAMS WEIBULL / Dongkook: You can ignore below - Just keep here
    PARAMS.CanStruc.beta(ii) = 6.6533;%2.1681;
    PARAMS.CanStruc.alpha(ii) = 0.8017;%0.5831;
    
    % Smax = maximum h2o storage capacity for foliage [mm/LAI unit]
    if isempty(cell2mat(para_leaf_crop_s(1,2)))
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop_s(1,3));
    else
        PARAMS.CanStruc.Smax = cell2mat(para_leaf_crop_s(1,2));
    end
    
    
    load (fullpath_P,...
        'para_canopy_crop_fixed')
    if isempty(cell2mat(para_canopy_crop_fixed(1,2)))
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,3));      % max fraction of canopy that can be wet
    else
        PARAMS.CanStruc.Ffact = cell2mat(para_canopy_crop_fixed(1,2));
    end
    if isempty(cell2mat(para_canopy_crop_fixed(2,2)))
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,3)); % precip extinction coefficient
    else
        PARAMS.CanStruc.pptintfact = cell2mat(para_canopy_crop_fixed(2,2));
    end
    
    % PHOTOSYNTHESIS
    
    load('./Temps/temp_model.mat', 'working_forcings', 'Start_Y')
    
    load (working_forcings)
    
    file_name_Py   = sprintf('YrwdSpec_%sY%s.mat', num2str(ii), num2str(Start_Y));
    
    load('.\Temps\temporary.mat', 'config_path');
    
    YrPath=fullfile(config_path, 'YearFr');
    
    fullpath_Py    =...
        fullfile(YrPath, file_name_Py);
    
    load (fullpath_Py, 'Ph_Type')
    
    load (fullpath_P,...
        'para_canopy_crop', 'para_photosynthesisC3_crop', 'para_photosynthesisC4_crop')
    %
    %     if ph_type==false
    %         disp('ok')
    %     end
    
    PARAMS.Photosyn.ph_type(ii) = Ph_Type;
    
    if PARAMS.Photosyn.ph_type(ii) == 1
        % C3 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC3_crop(1,2)))
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop(1,3));
        else
            PARAMS.Photosyn.beta_ph_C3(ii) = cell2mat(para_photosynthesisC3_crop(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop(2,2)))
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop(2,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Vcmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop(2,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop(3,2)))
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop(3,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Jmax25_C3{ii} = cell2mat(para_photosynthesisC3_crop(3,2))*ones(1,size(Ta_crop,1));
        end
        if isempty(cell2mat(para_photosynthesisC3_crop(4,2)))
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop(4,3))*ones(1,size(Ta_crop,1)); % [umol / m^2 / s]
        else
            PARAMS.Photosyn.Rd25{ii} = cell2mat(para_photosynthesisC3_crop(4,2))*ones(1,size(Ta_crop,1));
        end
        
        if isempty(cell2mat(para_photosynthesisC3_crop(5,2)))
            PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC3_crop(5,3)); % Vertical Distribution of Photosynthetic Capacity
        else
            PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC3_crop(5,2));
        end
        
        % C4 Photosynthesis Parameters
        PARAMS.Photosyn.Vmax_C4(ii) = NaN;
        PARAMS.Photosyn.Rd_C4(ii) = NaN;
        PARAMS.Photosyn.Q10_C4(ii) = NaN;
        PARAMS.Photosyn.kk_C4(ii) = NaN;
        PARAMS.Photosyn.theta_C4(ii) = NaN;
        PARAMS.Photosyn.beta_C4(ii) = NaN;
        PARAMS.Photosyn.al_C4(ii) = NaN;
    else
        % C3 Photosynthesis Parameters
        PARAMS.Photosyn.beta_ph_C3(ii) = NaN;
        PARAMS.Photosyn.Vcmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Jmax25_C3{ii} = NaN;
        PARAMS.Photosyn.Rd25{ii} = NaN;
        
        % C4 Photosynthesis Parameters
        if isempty(cell2mat(para_photosynthesisC4_crop(1,2)))
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop(1,3));
        else
            PARAMS.Photosyn.Vmax_C4(ii) = cell2mat(para_photosynthesisC4_crop(1,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop(2,2)))
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop(2,3));
        else
            PARAMS.Photosyn.Rd_C4(ii) = cell2mat(para_photosynthesisC4_crop(2,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop(3,2)))
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop(3,3));
        else
            PARAMS.Photosyn.Q10_C4(ii) = cell2mat(para_photosynthesisC4_crop(3,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop(4,2)))
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop(4,3));
        else
            PARAMS.Photosyn.kk_C4(ii) = cell2mat(para_photosynthesisC4_crop(4,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop(5,2)))
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop(5,3));
        else
            PARAMS.Photosyn.theta_C4(ii) = cell2mat(para_photosynthesisC4_crop(5,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop(6,2)))
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop(6,3));
        else
            PARAMS.Photosyn.beta_C4(ii) = cell2mat(para_photosynthesisC4_crop(6,2));
        end
        if isempty(cell2mat(para_photosynthesisC4_crop(7,2)))
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop(7,3));
        else
            PARAMS.Photosyn.al_C4(ii) = cell2mat(para_photosynthesisC4_crop(7,2));
        end
        
        if isempty(cell2mat(para_photosynthesisC4_crop(8,2)))
            PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC4_crop(8,3)); % Vertical Distribution of Photosynthetic Capacity
        else
            PARAMS.Photosyn.kn_canopy(1) = cell2mat(para_photosynthesisC4_crop(8,2));
        end
    end
    
    
    PARAMS.Photosyn.Oi(ii) = 210;                                           % Intercellular oxygen concentration [mmol / mol]
    
    
    % STOMATAL CONDUCTANCE
    load (fullpath_P,...
        'para_conductance_crop')
    
    % Ball-Berry
    if isempty(cell2mat(para_conductance_crop(1,2)))
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop(1,3)); % slope parameter in BB model [-] (Leakey: m=10.6 (ambient); m=10.9 (elevated))
    else
        PARAMS.StomCond.mslope(ii) = cell2mat(para_conductance_crop(1,2));
    end
    if isempty(cell2mat(para_conductance_crop(2,2)))
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop(2,3));   % intercept parameter in BB model [mol/m^2/s] (Leakey: b=0.008 (ambient); b=0.007 (elevated))
    else
        PARAMS.StomCond.bint(ii) = cell2mat(para_conductance_crop(2,2));
    end
    
    
    % (Tuzet et al, PCE 2003)
    if isempty(cell2mat(para_conductance_crop(3,2)))
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop(3,3));    % sensitivity parameter for initial decrease in leaf potential [-]
    else
        PARAMS.StomCond.sfm(ii) = cell2mat(para_conductance_crop(3,2));
    end
    if isempty(cell2mat(para_conductance_crop(4,2)))
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop(4,3));  % leaf potential at which half of the hydraulic conductance is lost [MPa]
    else
        PARAMS.StomCond.psifm(ii) = cell2mat(para_conductance_crop(4,2));
    end
    
    
    % Conductivities
    if isempty(cell2mat(para_conductance_crop(5,2)))
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop(5,3));      % radial conductivity of the root system [s^-1]
    else
        PARAMS.Soil.K_rad(ii) = cell2mat(para_conductance_crop(5,2));
    end
    if isempty(cell2mat(para_conductance_crop(6,2)))
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop(6,3));      % axial specific conductivity of the root system [mm/s]
    else
        PARAMS.Soil.K_axs(ii) = cell2mat(para_conductance_crop(6,2));
    end
    if isempty(cell2mat(para_conductance_crop(7,2)))
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop(7,3));
    else
        PARAMS.Soil.kpar_ax(ii) = cell2mat(para_conductance_crop(7,2));
    end
    
    if isempty(cell2mat(para_conductance_crop(8,2)))
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop(8,3));
    else
        PARAMS.StomCond.Rp(ii) = cell2mat(para_conductance_crop(8,2));
    end
    
    % ROOT SOIL
    % HR
    
    load('./Temps/temp_model.mat', 'HR')
    
    SWITCHES.HR_on(ii) = HR;                                          % Hydraulic Redistribution 1. Yes, 0. No
    
    % root structure
    load (fullpath_Py,...
        'set_para_root')
    PARAMS.Soil.z50f(ii) = cell2mat(set_para_root(2,2));             % From Ameriflux site description
    PARAMS.Soil.z95f(ii) = cell2mat(set_para_root(3,2));             % From Ameriflux site description
    PARAMS.Soil.z50t(ii) = cell2mat(set_para_root(2,2));
    PARAMS.Soil.z95t(ii) = cell2mat(set_para_root(3,2));
    PARAMS.Soil.maxrootdepth(ii) = cell2mat(set_para_root(1,2));
    
    
    % Cropgrowth
    
    load(fullpath_P, 'para_CGMF')
    
    % Fraction of carbohydrate flux to leaf
    if isempty(cell2mat(para_CGMF(1,2)))
        PARAMS.CGM.LFPT_3(ii)=cell2mat(para_CGMF(1,3));
    else
        PARAMS.CGM.LFPT_3(ii)=cell2mat(para_CGMF(1,2));
    end
    
    if isempty(cell2mat(para_CGMF(2,2)))
        PARAMS.CGM.LFPT_4(ii)=cell2mat(para_CGMF(2,3));
    else
        PARAMS.CGM.LFPT_4(ii)=cell2mat(para_CGMF(2,2));
    end
    
    if isempty(cell2mat(para_CGMF(3,2)))
        PARAMS.CGM.LFPT_5(ii)=cell2mat(para_CGMF(3,3));
    else
        PARAMS.CGM.LFPT_5(ii)=cell2mat(para_CGMF(3,2));
    end
    
    if isempty(cell2mat(para_CGMF(4,2)))
        PARAMS.CGM.LFPT_6(ii)=cell2mat(para_CGMF(4,3));
    else
        PARAMS.CGM.LFPT_6(ii)=cell2mat(para_CGMF(4,2));
    end
    
    
    % Fraction of carbohydrate flux to stem
    if isempty(cell2mat(para_CGMF(5,2)))
        PARAMS.CGM.STPT_3(ii)=cell2mat(para_CGMF(5,3));
    else
        PARAMS.CGM.STPT_3(ii)=cell2mat(para_CGMF(5,2));
    end
    
    if isempty(cell2mat(para_CGMF(6,2)))
        PARAMS.CGM.STPT_4(ii)=cell2mat(para_CGMF(6,3));
    else
        PARAMS.CGM.STPT_4(ii)=cell2mat(para_CGMF(6,2));
    end
    
    if isempty(cell2mat(para_CGMF(7,2)))
        PARAMS.CGM.STPT_5(ii)=cell2mat(para_CGMF(7,3));
    else
        PARAMS.CGM.STPT_5(ii)=cell2mat(para_CGMF(7,2));
    end
    
    if isempty(cell2mat(para_CGMF(8,2)))
        PARAMS.CGM.STPT_6(ii)=cell2mat(para_CGMF(8,3));
    else
        PARAMS.CGM.STPT_6(ii)=cell2mat(para_CGMF(8,2));
    end
    
    
    % Fraction of carbohydrate flux to root
    
    if isempty(cell2mat(para_CGMF(9,2)))
        PARAMS.CGM.RTPT_3(ii)=cell2mat(para_CGMF(9,3));
    else
        PARAMS.CGM.RTPT_3(ii)=cell2mat(para_CGMF(9,2));
    end
    
    if isempty(cell2mat(para_CGMF(10,2)))
        PARAMS.CGM.RTPT_4(ii)=cell2mat(para_CGMF(10,3));
    else
        PARAMS.CGM.RTPT_4(ii)=cell2mat(para_CGMF(10,2));
    end
    
    if isempty(cell2mat(para_CGMF(11,2)))
        PARAMS.CGM.RTPT_5(ii)=cell2mat(para_CGMF(11,3));
    else
        PARAMS.CGM.RTPT_5(ii)=cell2mat(para_CGMF(11,2));
    end
    
    if isempty(cell2mat(para_CGMF(12,2)))
        PARAMS.CGM.RTPT_6(ii)=cell2mat(para_CGMF(12,3));
    else
        PARAMS.CGM.RTPT_6(ii)=cell2mat(para_CGMF(12,2));
    end
    
    % Fraction of carbohydrate flux to grain
    
    if isempty(cell2mat(para_CGMF(13,2)))
        PARAMS.CGM.GRAINPT_5(ii)=cell2mat(para_CGMF(13,3));
    else
        PARAMS.CGM.GRAINPT_5(ii)=cell2mat(para_CGMF(13,2));
    end
    
    if isempty(cell2mat(para_CGMF(14,2)))
        PARAMS.CGM.GRAINPT_6(ii)=cell2mat(para_CGMF(14,3));
    else
        PARAMS.CGM.GRAINPT_6(ii)=cell2mat(para_CGMF(14,2));
    end
    
    
    load(fullpath_P, 'para_CGM1')
    
    % Base temperature for GDD accumulation
    if isempty(cell2mat(para_CGM1(1,2)))
        PARAMS.CGM.GDDTBASE(ii)=cell2mat(para_CGM1(1,3));
    else
        PARAMS.CGM.GDDTBASE(ii)=cell2mat(para_CGM1(1,2));
    end
    
    % Upper temperature for GDD accumulation (c)
    if isempty(cell2mat(para_CGM1(2,2)))
        PARAMS.CGM.GDDTCUT(ii)=cell2mat(para_CGM1(2,3));
    else
        PARAMS.CGM.GDDTCUT(ii)=cell2mat(para_CGM1(2,2));
    end
    
    % GDD from seeding to emergence (upto VE)
    if isempty(cell2mat(para_CGM1(3,2)))
        PARAMS.CGM.GDDS1(ii)=cell2mat(para_CGM1(3,3));
    else
        PARAMS.CGM.GDDS1(ii)=cell2mat(para_CGM1(3,2));
    end
    
    % GDD from seeding to initial vegetative (upto V15)
    if isempty(cell2mat(para_CGM1(4,2)))
        PARAMS.CGM.GDDS2(ii)=cell2mat(para_CGM1(4,3));
    else
        PARAMS.CGM.GDDS2(ii)=cell2mat(para_CGM1(4,2));
    end
    
    % GDD from seeding to normal vegetative (upto VT)
    if isempty(cell2mat(para_CGM1(5,2)))
        PARAMS.CGM.GDDS3(ii)=cell2mat(para_CGM1(5,3));
    else
        PARAMS.CGM.GDDS3(ii)=cell2mat(para_CGM1(5,2));
    end
    
    % GDD from seeding to initial reproductive (upto R2)
    if isempty(cell2mat(para_CGM1(6,2)))
        PARAMS.CGM.GDDS4(ii)=cell2mat(para_CGM1(6,3));
    else
        PARAMS.CGM.GDDS4(ii)=cell2mat(para_CGM1(6,2));
    end
    
    % GDD from seeding to physical maturity (upto R6)
    if isempty(cell2mat(para_CGM1(7,2)))
        PARAMS.CGM.GDDS5(ii)=cell2mat(para_CGM1(7,3));
    else
        PARAMS.CGM.GDDS5(ii)=cell2mat(para_CGM1(7,2));
    end
    
    % Q10 for maintenance respiration
    if isempty(cell2mat(para_CGM1(8,2)))
        PARAMS.CGM.Q10MR(ii)=cell2mat(para_CGM1(8,3));
    else
        PARAMS.CGM.Q10MR(ii)=cell2mat(para_CGM1(8,2));
    end
    
    
    % Leaf area per living leaf biomass (SLA) By 2013 ISAM papers
    if isempty(cell2mat(para_CGM1(9,2)))
        PARAMS.CGM.BIO2LAI(ii)=cell2mat(para_CGM1(9,3));
    else
        PARAMS.CGM.BIO2LAI(ii)=cell2mat(para_CGM1(9,2));
    end
    
    % Leaf maintenance respiration at 25c (µmol co2/m**2/s)
    if isempty(cell2mat(para_CGM1(10,2)))
        PARAMS.CGM.LFMR25(ii)=cell2mat(para_CGM1(10,3));
    else
        PARAMS.CGM.LFMR25(ii)=cell2mat(para_CGM1(10,2));
    end
    
    % Stem maintenance respiration at 25c (µmol co2/kg bio/s)
    if isempty(cell2mat(para_CGM1(11,2)))
        PARAMS.CGM.STMR25(ii)=cell2mat(para_CGM1(11,3));
    else
        PARAMS.CGM.STMR25(ii)=cell2mat(para_CGM1(11,2));
    end
    
    % Root maintenance respiration at 25c (µmol co2/kg bio/s)
    if isempty(cell2mat(para_CGM1(12,2)))
        PARAMS.CGM.RTMR25(ii)=cell2mat(para_CGM1(12,3));
    else
        PARAMS.CGM.RTMR25(ii)=cell2mat(para_CGM1(12,2));
    end
    
    % Grain maintenance respiration at 25c (µmol co2/kg bio/s)
    if isempty(cell2mat(para_CGM1(13,2)))
        PARAMS.CGM.GRAINMR25(ii)=cell2mat(para_CGM1(13,3));
    else
        PARAMS.CGM.GRAINMR25(ii)=cell2mat(para_CGM1(13,2));
    end
    
    % Fraction of growth respiration
    if isempty(cell2mat(para_CGM1(14,2)))
        PARAMS.CGM.FRA_GR(ii)=cell2mat(para_CGM1(14,3));
    else
        PARAMS.CGM.FRA_GR(ii)=cell2mat(para_CGM1(14,2));
    end
    
    % Leaf turnover coefficient (1/s)
    if isempty(cell2mat(para_CGM1(15,2)))
        PARAMS.CGM.LF_OVRC_5(ii)=cell2mat(para_CGM1(15,3));
    else
        PARAMS.CGM.LF_OVRC_5(ii)=cell2mat(para_CGM1(15,2));
    end
    
    if isempty(cell2mat(para_CGM1(16,2)))
        PARAMS.CGM.LF_OVRC_6(ii)=cell2mat(para_CGM1(16,3));
    else
        PARAMS.CGM.LF_OVRC_6(ii)=cell2mat(para_CGM1(16,2));
    end
    
    % Stem turnover coefficient (1/s)
    if isempty(cell2mat(para_CGM1(17,2)))
        PARAMS.CGM.ST_OVRC_5(ii)=cell2mat(para_CGM1(17,3));
    else
        PARAMS.CGM.ST_OVRC_5(ii)=cell2mat(para_CGM1(17,2));
    end
    
    if isempty(cell2mat(para_CGM1(18,2)))
        PARAMS.CGM.ST_OVRC_6(ii)=cell2mat(para_CGM1(18,3));
    else
        PARAMS.CGM.ST_OVRC_6(ii)=cell2mat(para_CGM1(18,2));
    end
    
    % Root turnover coefficient (1/s)
    if isempty(cell2mat(para_CGM1(19,2)))
        PARAMS.CGM.RT_OVRC_5(ii)=cell2mat(para_CGM1(19,3));
    else
        PARAMS.CGM.RT_OVRC_5(ii)=cell2mat(para_CGM1(19,2));
    end
    
    if isempty(cell2mat(para_CGM1(20,2)))
        PARAMS.CGM.RT_OVRC_6(ii)=cell2mat(para_CGM1(20,3));
    else
        PARAMS.CGM.RT_OVRC_6(ii)=cell2mat(para_CGM1(20,2));
    end
    
    if isempty(cell2mat(para_CGM1(21,2)))
        PARAMS.CGM.CF_LF(ii)=cell2mat(para_CGM1(21,3));
    else
        PARAMS.CGM.CF_LF(ii)=cell2mat(para_CGM1(21,2));
    end
    
    if isempty(cell2mat(para_CGM1(22,2)))
        PARAMS.CGM.CF_ST(ii)=cell2mat(para_CGM1(22,3));
    else
        PARAMS.CGM.CF_ST(ii)=cell2mat(para_CGM1(22,2));
    end
    
    if isempty(cell2mat(para_CGM1(23,2)))
        PARAMS.CGM.CF_RT(ii)=cell2mat(para_CGM1(23,3));
    else
        PARAMS.CGM.CF_RT(ii)=cell2mat(para_CGM1(23,2));
    end
    
    if isempty(cell2mat(para_CGM1(24,2)))
        PARAMS.CGM.CF_GN(ii)=cell2mat(para_CGM1(24,3));
    else
        PARAMS.CGM.CF_GN(ii)=cell2mat(para_CGM1(24,2));
    end
    
    if isempty(cell2mat(para_CGM1(25,2)))
        PARAMS.CGM.HI(ii)=cell2mat(para_CGM1(25,3));
    else
        PARAMS.CGM.HI(ii)=cell2mat(para_CGM1(25,2));
    end
    
    % Foliage nitrogen concentration (%)
    if isempty(cell2mat(para_CGM1(26,2)))
        PARAMS.CGM.FOLN_MX(ii)=cell2mat(para_CGM1(26,3));
    else
        PARAMS.CGM.FOLN_MX(ii)=cell2mat(para_CGM1(26,2));
    end
    
    % Characteristic T for leaf freezing (K)
    if isempty(cell2mat(para_CGM1(27,2)))
        PARAMS.CGM.LEFREEZ(ii)=cell2mat(para_CGM1(27,3));
    else
        PARAMS.CGM.LEFREEZ(ii)=cell2mat(para_CGM1(27,2));
    end
    
    % Coefficient for leaf temperature stress death (1/s)
    if isempty(cell2mat(para_CGM1(28,2)))
        PARAMS.CGM.DILE_FC_5(ii)=cell2mat(para_CGM1(28,3));
    else
        PARAMS.CGM.DILE_FC_5(ii)=cell2mat(para_CGM1(28,2));
    end
    
    if isempty(cell2mat(para_CGM1(29,2)))
        PARAMS.CGM.DILE_FC_6(ii)=cell2mat(para_CGM1(29,3));
    else
        PARAMS.CGM.DILE_FC_6(ii)=cell2mat(para_CGM1(29,2));
    end
    
    
    PARAMS.CGM.DILE_FW_5 = 0;                                 % Coefficient for leaf water stress death (1/s)
    PARAMS.CGM.DILE_FW_6 = 0;
    PARAMS.CGM.WSTRES = 0;                                                 % Soil water stress
    
    
    % N-allocation
    
    load(fullpath_P, 'para_N_allocation')
    
    
    if isempty(cell2mat(para_N_allocation(1,2)))
        PARAMS.Nalloc.Nstr(ii)=cell2mat(para_N_allocation(1,3));
    else
        PARAMS.Nalloc.Nstr(ii)=cell2mat(para_N_allocation(1,2));
    end
    
    if isempty(cell2mat(para_N_allocation(2,2)))
        PARAMS.Nalloc.Nlhc(ii)=cell2mat(para_N_allocation(2,3));
    else
        PARAMS.Nalloc.Nlhc(ii)=cell2mat(para_N_allocation(2,2));
    end
    
    if isempty(cell2mat(para_N_allocation(3,2)))
        PARAMS.Nalloc.Rnchl(ii)=cell2mat(para_N_allocation(3,3));
    else
        PARAMS.Nalloc.Rnchl(ii)=cell2mat(para_N_allocation(3,2));
    end
    
    if isempty(cell2mat(para_N_allocation(4,2)))
        PARAMS.Nalloc.MinKec(ii)=cell2mat(para_N_allocation(4,3));
    else
        PARAMS.Nalloc.MinKec(ii)=cell2mat(para_N_allocation(4,2));
    end
    
    if isempty(cell2mat(para_N_allocation(5,2)))
        PARAMS.Nalloc.MaxKec(ii)=cell2mat(para_N_allocation(5,3));
    else
        PARAMS.Nalloc.MaxKec(ii)=cell2mat(para_N_allocation(5,2));
    end
    
    if isempty(cell2mat(para_N_allocation(6,2)))
        PARAMS.Nalloc.MinFth(ii)=cell2mat(para_N_allocation(6,3));
    else
        PARAMS.Nalloc.MinFth(ii)=cell2mat(para_N_allocation(6,2));
    end
    
    if isempty(cell2mat(para_N_allocation(7,2)))
        PARAMS.Nalloc.MaxFth(ii)=cell2mat(para_N_allocation(7,3));
    else
        PARAMS.Nalloc.MaxFth(ii)=cell2mat(para_N_allocation(7,2));
    end
    
    if isempty(cell2mat(para_N_allocation(8,2)))
        PARAMS.Nalloc.MinFrub(ii)=cell2mat(para_N_allocation(8,3));
    else
        PARAMS.Nalloc.MinFrub(ii)=cell2mat(para_N_allocation(8,2));
    end
    
    if isempty(cell2mat(para_N_allocation(9,2)))
        PARAMS.Nalloc.MaxFrub(ii)=cell2mat(para_N_allocation(9,3));
    else
        PARAMS.Nalloc.MaxFrub(ii)=cell2mat(para_N_allocation(9,2));
    end
    
end




% Layer to cut
load (fullpath_Py,...
    'set_para_root')
CONSTANTS.nlc = set_para_root{4,2};            % Number of layer to cut
VARIABLES.comroot = 0;                          % Initialize to zero the cut of roots


if set_para_root{4,2} > 0
    SWITCHES.cutroots = 1;        %  0 if no cut at all
    %  (nlc) 1 if cut all the time same layers
    %  2 if cut with a threshold
    %  3 if cut with a embolism curve (PLC)
    %  4 combination of 1, 2 and 3
elseif set_para_root{4,2} == 0
    SWITCHES.cutroots = 0;        %  0 if no cut at all
    %  (nlc) 1 if cut all the time same layers
    %  2 if cut with a threshold
    %  3 if cut with a embolism curve (PLC)
    %  4 combination of 1, 2 and 3
end


%*************************************************************************%
%********************** Site Specific Parameters *************************%
%*************************************************************************%
load (fullpath_Py,...
    'num_LAD')
if isstr(num_LAD) == 1
    PARAMS.CanStruc.nl_can=str2num(num_LAD);   % # canopy layers
else
    PARAMS.CanStruc.nl_can=(num_LAD);          % # canopy layers
end


%*************************************************************************%
%**************************** Flux / Met Data ****************************%
%******************************************************* *****************%

% LOAD DATA
decyear = nan;

load './Temps/temp_model.mat'...
    'working_forcings' 'lat_face' 'long_face' 'elev_face'
load (working_forcings)


% Organize the data
% 
% if length(doys)==365
%     inds = find(ismember(year_crop, Run_years));
% else
%     inds = find(ismember(year_crop, Run_years) & ismember(floor(round((doy_crop+1).*10^10)./10^10), doys));
% end
    


inds = find(ismember(year_crop, Run_years) & ismember(floor(round((doy_crop+1).*10^10)./10^10), doys));

ELEV      = elev_face;
LAT       = lat_face;
LONG      = long_face;
decyear   = year_crop(inds);
decdoy    = doy_crop(inds);

year      = year_crop(inds);
doy       = doy_crop(inds);
hour      = hour_crop(inds);
ZEN_in    = ZEN_crop(inds);
Rg_in     = Rg_crop(inds);

Ta_in     = Ta_crop(inds);
VPD_in    = VPD_crop(inds);
PPT_in    = PPT_crop(inds);
U_in      = U_crop(inds);
ustar_in  = ustar_crop(inds);

Pa_in     = Pa_crop(inds);
ea_in     = ea_crop(inds);

Ta_all     = Ta_crop;
management = management_crop;

% Dongkook: Calucalated in model_forcing
% % Calculate Vapor Variables
%     aa = 0.611;  % [kPa]
%     bb = 17.502;
%     cc = 240.97; % [C]
%     esat = aa*exp((bb*Ta)./(cc+Ta));
%     ean = esat - VPD;
%     hr = ea./esat;
%

% Data Corrections
uinds = find(U_in<.1); %1
U_in(uinds) = .1;

% Calculate ustar for missing periods
binds = find(isnan(ustar_in) | ustar_in<=0);
vonk = CONSTANTS.vonk;
hobs = PARAMS.CanStruc.hobs;
z0 = PARAMS.CanStruc.z0;
hcan =PARAMS.CanStruc.hcan;
d0 = PARAMS.CanStruc.d0;
ustar_in(binds) = vonk.*U_in(binds)./(log((hobs-d0)/z0));

% Adjust U from measurement height to canopy height (Brutsaert, 4.2)
U_in = U_in - (ustar_in/vonk).*log(hobs/hcan);
uinds = find(U_in<.1);
U_in(uinds) = .1;

% Dongkook: No idea where it is used
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% % FILL Measured Hg with Linear fit of Hg to Rg
% % --> Using Hg to force soil heat transfer model
%     ginds = find(~isnan(Hg) & ~isnan(Rg));
%     [lps] = polyfit(Rg(ginds), Hg(ginds), 1);
%     aa = lps(1); bb = lps(2);
%     binds = find(isnan(Hg));
%     Hg(binds) = aa*Rg(binds) + bb;
% %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


if (kspecies == 1)
    LAI_in(:,1)= LAI_species1_crop(inds);
elseif (kspecies == 2)
    LAI_in(:,1)= LAI_species1_crop(inds);
    LAI_in(:,2)= LAI_species2_crop(inds);
elseif (kspecies == 3)
    LAI_in(:,1)= LAI_species1_crop(inds);
    LAI_in(:,2)= LAI_species2_crop(inds);
    LAI_in(:,3)= LAI_species3_crop(inds);
elseif (kspecies == 4)
    LAI_in(:,1)= LAI_species1_crop(inds);
    LAI_in(:,2)= LAI_species2_crop(inds);
    LAI_in(:,3)= LAI_species3_crop(inds);
    LAI_in(:,4)= LAI_species4_crop(inds);
elseif (kspecies == 5 && Sim_species_con == 1)
    LAI_in(:,1)= LAI_species1_crop(inds);
elseif (kspecies == 5 && Sim_species_con == 2)
    LAI_in(:,1)= LAI_species2_crop(inds);
elseif (kspecies == 5 && Sim_species_con == 3)
    LAI_in(:,1)= LAI_species3_crop(inds);
elseif (kspecies == 5 && Sim_species_con == 4)
    LAI_in(:,1)= LAI_species4_crop(inds);
end

nspecies=PARAMS.CanStruc.nspecies;

% ALLOCATE MEMEORY FOR NVINDS
nvinds_all = cell(nspecies,1);
vinds_all = cell(nspecies,1);

