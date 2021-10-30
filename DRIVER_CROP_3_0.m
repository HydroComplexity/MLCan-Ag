function DRIVER_CROP_3_0
% Canopy-Root-Soil-Atmosphere Exchange Model
%
%   Written By : Darren Drewry (dtd2@illinois.edu)
%              : Juan Quijano (quijano2@illinois.edu)
%              : Dongkook Woo (dwoo5@illinois.edu)

% STRUCTURES:
%   SWITCHES --> model conditional switches
%   VERTSTRUC --> variables describing vertical structure of canopy & soil
%   VARIABLES.CANOPY --> canopy variables
%            .SOIL --> soil variables
%            .ROOT --> root variables
%   FORCING --> holds current timestep forcing variables
%   CONSTANTS --> site independent constants, unit conversions
%   PARAMS
%         .CanStruc --> canopy structural parameters
%         .Rad --> radiation parameters
%         .Photosyn --> photosynthesis paramters
%         .StomCond --> stomtatal conductance parameters
%         .Resp --> ecosystem respiration parameters
%         .MicroEnv --> canopy microenvironment parameters
%         .Soil --> soil paramters
%
%**************************************************************************
%                           USER SPECIFICATIONS
%**************************************************************************
load './Temps/temp_model.mat' ...
    'DOY_start' 'DOY_end'

doys = [DOY_start:DOY_end];

CROP_GROWTH.DOY_start = DOY_start; %Rohit added
CROP_GROWTH.DOY_end = DOY_end; %Rohit added

%******************************  SWITCHES   *******************************
load ('./Temps/temp_model.mat',...
    'CGM', 'Turbulence', 'HR', 'Soil_heat', 'Soil_nutrient', 'NCanAlMod', 'ParOn',...
    ...
    'vanGen', 'RHC', 'CO2_Ambient',...
    'CO2_Elev', 'CO2_Elev_con', 'Temp_Elev', 'Temp_Elev_con')



SWITCHES.plants = 0;              %  Allow the presence or no of plants. 1. With plants 2. No plants (bare soil, no roots)
SWITCHES.LT = 0;                  %  Run long term dynamics with stochastically generated data 1. Yes 0. No
SWITCHES.turb_on = Turbulence;    %  1 = scalar profiles resolved, otherwise not resolved

SWITCHES.littertype = 2;          %  [1 LEs parameterized with aerodynamic resistance
%  2 LEs parameterized with vapor litter diffusivity]

SWITCHES.soilevap = 1;            %  1 = Include soil evaporation. 0 does not include soil evaporation

SWITCHES.rtcond  = 0;             %  0. Juan's Method, 1. Amenu Method

SWITCHES.vanGen = vanGen;         %  1. vanGenuchten, 0. Brooks Corey / Soil moisture chracteristic.

SWITCHES.Pedofunctions = 0;       %  1. Use pedo transfer functions with organic matter 0. Do not use

SWITCHES.rhc = RHC;               %  1 = linearly increasing root hydraulic conductivity with depth

SWITCHES.ns=0;                    %  Numerical Scheme. 1 = Use Implicit, 0 = Use explicit

SWITCHES.soilheat_on = Soil_heat; %  Compute the heat equation 1 = Yes, 0 = No

SWITCHES.canstorheat = 1;         %  Include storage of heat in leaves 1 = Yes, 0 = No
SWITCHES.useG_on=0;               %  Use of G measured in the soil energy balance instead of compute it 0 = Off, 1 = On
SWITCHES.useTs_on=0;              %  Use of Ts (soil temperature) [soil energy balance]

SWITCHES.entropy_on = 0;          %  Compute Entropy, 1 = Yes, 0 = No.
SWITCHES.entropymethod = 2;       %  Compute Entropy Method for LW entropy. 1. Landsberg 1979. 2. Smith 2001. 3. Clausius ***MEREDITH: 3 nonexistent in code; 7/26/2017***

SWITCHES.fsv_off = 0;             %  1 = hydraulic constraint turned OFF

SWITCHES.temp_change = Temp_Elev; %  0 = +0 C, 1 = 1 C,  2 = 2 C

% Added by RN
SWITCHES.Onset = 0;     % For CGM 1- Select planting date based on soil moisture
% 0- use the provided data of planting date
SWITCHES.CGM = CGM;
SWITCHES.NCanAlMod = NCanAlMod;

if strcmpi(ParOn, 'On')
    SWITCHES.ParOnV = true;
else
    SWITCHES.ParOnV = false;
end

% All of the CN model switch and parameters are in core_N
SWITCHES.soilCN_on = Soil_nutrient;%  Compute CN dynamics. 1 = Yes, 0 = No.


SWITCHES.save_on = 1;             %  1 = save stored variables to .mat file, otherwise no save performed
SWITCHES.plots_on = 0;

%**************************************************************************
% Code Library Paths
addpath('./LOCAL_CODES/CANOPY/');
addpath('./LOCAL_CODES/ENTROPY/');
addpath('./LOCAL_CODES/ROOT_SOIL/');
addpath('./LOCAL_CODES/ROOT_SOIL/IMPLICIT');
addpath('./LOCAL_CODES/ROOT_SOIL/CN_MODEL');
addpath('./LOCAL_CODES/NUMERICAL/');
addpath('./LOCAL_CODES/NUMERICAL/OTHERS');
addpath('./LOCAL_CODES/PLOTTING/');
addpath('./LOCAL_CODES/CROP_GROWTH/');
addpath('./LOCAL_CODES/IRRIGATION/');
addpath('./LOCAL_CODES/CANOPY_NITROGEN/');



load('./Temps/temp_model.mat', 'Start_Y', 'End_Y')

Run_years=Start_Y:End_Y;

%************************** LOAD INFORMATION  *****************************
year = nan;                       % Initialize year, otherwise it is recognized as function
LOAD_SITE_INFO;


% ybeginds = nan(length(Each_year),1);
% yendinds = nan(length(Each_year),1);


for yy = 1:length(Run_years)
    ybeginds(yy) = find(year==Run_years(yy), 1, 'first');
    yendinds(yy) = find(year==Run_years(yy), 1, 'last');
end

%************************** CLIMATE CHANGE ********************************
if (SWITCHES.temp_change == 0)    % Ambient Temperature
    temp_change = 0;
elseif SWITCHES.temp_change == 1  % Elevated Temperature
    temp_change = Temp_Elev_con;
    Ta_in = Ta_in + temp_change;
end

CO2base_elevated = CO2_Elev_con;

CO2base_ambient = CO2_Ambient;

if (CO2_Elev == 1)                % Elevated CO2
    CO2base = CO2base_elevated;
else                              % Ambient CO2
    CO2base = CO2base_ambient;
end

%************************ CANOPY - SOIL - FUNCTION ************************
CANOPY_SOIL_COUPLER;
%**************************************************************************
%
% workpath=pwd;
%
% cd(workpath);

msgbox(['Simulation is done.', ResFile, ' is saved in the Results folder'],'MLCan Simulation');


%**************************************************************************
