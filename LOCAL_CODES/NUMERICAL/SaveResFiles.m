
%******************************* SAVE FUNCTION ****************************
load './Temps/temp_model.mat'...
    'num_species'
if isstr(num_species) == 1
    Number_S = num_species;
else
    Number_S = num2str(num_species);
end
current_time=clock;
current_year=num2str(current_time(1));
current_month=num2str(current_time(2));
current_day=num2str(current_time(3));
current_hour=num2str(current_time(4));
current_min=num2str(current_time(5));
simulation_year=num2str(Run_years(yy));

% Results_folder_name = [];
%
% load './Temps/temp_model.mat' ...
%     'Results_folder_name' % Rohit added

%% Rohit Start
%
% if  isempty(Results_folder_name)
%     Results_folder_name = '.\Result\';
% end

load('.\Temps\temporary.mat', 'config_path');

ResPath=fullfile(config_path, 'Results');

Results_folder_name=ResPath;

ResFile=['Result_', Number_S, 'species', '_sim_year_', simulation_year, '_Saved_on_' ...
    , current_month, '.', current_day, '.', current_year, '.mat'];

namefileout = fullfile(Results_folder_name, ResFile);

% Make dianuual averages
MAKE_AVERAGES;

% % Show message
% if (length(Each_year)==1)
%     timevect  = decdoy;
%     timelabel = ['DOY: ', num2str(Each_year)];
% else
%     timevect  = [1:length(decyear)];
%     %       timevect  = decyear;
%     timelabel = ['Decimal Year'];
% end

plantingDate= CROP_GROWTH.plantingDate;
doy1 = CROP_GROWTH.doy1;
doy2 = CROP_GROWTH.doy2;
doy3 = CROP_GROWTH.doy3;
doy4 = CROP_GROWTH.doy4;
doy5 = CROP_GROWTH.doy5;
harvestingDate= CROP_GROWTH.harvestingDate;


MLCan_selected_namefileout = fullfile(Results_folder_name,...
    ['Result_selected', '_sim_year_', simulation_year, '.mat']);
Biomassnamefileout = fullfile(Results_folder_name,...
    ['Result_BM', '_sim_year_', simulation_year, '.mat']);


save (namefileout)


if SWITCHES.soilCN_on
    save(MLCan_selected_namefileout, 'volliq_store', 'Tl_mean', 'Ts_store',...
        'An_canopy_prof', 'An_can_store', 'TR_can_all_store',...
        'H_canopy_prof', 'H_can_store', 'Evap_canopy', 'Evap_canopy_prof',...
        'Esoil_store', 'doy1', 'doy2', 'doy3', 'doy4', 'doy5',...
        'plantingDate', 'harvestingDate', 'Cl_store', 'Ch_store',...
        'Cb_store', 'Amm_store', 'Nit_store', 'CNl_store', 'UP_N_m2_store',...
        'Chl_la_store', 'Vcmax_vert_store', 'N_can_store', 'N_rub_la_store',...
        'N_leaf_la_store', 'N_uptake_can_store', 'K_ec_store',...
        'K_rub_leaf_store', 'K_lhc_leaf_store', 'psil_mean', 'fn_store', 'Nreq_store',...
        'Fert_req_store', 'Fert_doy_store')
    
else
    save(MLCan_selected_namefileout, 'volliq_store', 'Tl_mean', 'Ts_store', 'An_canopy_prof', 'An_can_store',...
        'TR_can_all_store', 'H_canopy_prof', 'H_can_store', 'Evap_canopy', 'Evap_canopy_prof', 'Esoil_store',...
        'doy1', 'doy2', 'doy3', 'doy4', 'doy5', 'plantingDate', 'harvestingDate')
end

if SWITCHES.CGM==1
    save(Biomassnamefileout, 'LeafCMass_can_store', 'StemCMass_can_store', 'RootCMass_can_store',...
        'GrainCMass_can_store', 'LAIsim_can_store', 'TotalCMass_can_store', 'AboveCMass_can_store',...
        'LeafBiomass_can_store', 'StemBiomass_can_store', 'RootBiomass_can_store', 'GrainBiomass_can_store',...
        'TotalBioMass_can_store', 'AboveBioMass_can_store', 'GrainYield_store', 'IW_store',...
        'NDemand_store', 'DeathRate_store')
end

%% Rohit End