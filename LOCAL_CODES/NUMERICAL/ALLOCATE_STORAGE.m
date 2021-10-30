% Script to allocate storage vectors / matrices

% N = length(doy);

N = NoDay*24*60/CONSTANTS.timestep;

% CANOPY VARIABLES
PARabs_sun_prof = NaN(nl_can, N);           
PARabs_shade_prof = NaN(nl_can, N);         
PARabs_canopy_prof = NaN(nl_can, N);        

PARabs_sun_norm_prof = NaN(nl_can, N);      
PARabs_shade_norm_prof = NaN(nl_can, N);    
PARabs_canopy_norm_prof = NaN(nl_can, N);   

NIRabs_sun_prof = NaN(nl_can, N);           
NIRabs_shade_prof = NaN(nl_can, N);         
NIRabs_canopy_prof = NaN(nl_can, N);        

NIRabs_sun_norm_prof = NaN(nl_can, N);      
NIRabs_shade_norm_prof = NaN(nl_can, N);    
NIRabs_canopy_norm_prof = NaN(nl_can, N);   

LWabs_can_prof = NaN(nl_can, N);            
LWemit_can_prof = NaN(nl_can, N);           

LWemit_can_store = NaN(1, N);           
LWemit_soil_store = NaN(1, N);          

SWout_store = NaN(1, N);           
LWout_store = NaN(1, N);           
LWdn_in = NaN(1, N);          
fdiff_store = NaN(1, N);      
refl_soil_store = NaN(1, N);  

% FLUXES:
% Ecosystem Fluxes (Canopy + Soil)
Fc_eco_store = NaN(1, N);     
LE_eco_store = NaN(1, N);     
H_eco_store = NaN(1, N);      
Rnrad_eco_store = NaN(1, N);            % Rnrad_eco;                                         %[W/m^2] Ecosystem Net Radiation

% Canopy Fluxes
Ph_can_store = NaN(1, N);            % Ph_can;                                               %[umol CO2/ m^2 ground / s] Photosynthetic flux from ecosystem
Ph_can_all_store = NaN(1, N, nspecies);                      %[umol CO2/ m^2 ground / s] Photosynthetic flux for each species
An_can_store = NaN(1, N);            % An_can;                                               %[umol CO2/ m^2 ground / s] Photosynthetic flux minus leaf respiration from ecosystem
An_can_all_store = NaN(1, N, nspecies);                      %[umol CO2/ m^2 ground / s] Photosynthetic flux minus from each species
LE_can_store = NaN(1, N);            % LE_can;                                               %[W/m^2] Total Latent Heat from Vegetation
LE_can_all_store = NaN(1, N, nspecies);                      %[W/m^2] Total Latent Heat from Vegetation for each species
H_can_store = NaN(1, N);            % H_can;                                                 %[W/m^2] Total Sensible Heat from Vegetation
H_can_all_store = NaN(1, N, nspecies);                        %[W/m^2] Total Sensible Heat from Vegetation for each species
dH_can_store = NaN(1, N);            % dHcan;                                                %[W/m^2] Total change of heat content in leaves
dH_can_all_store = NaN(1, N, nspecies);                      %[W/m^2] Total change of heat content in leaves for each species
TR_can_store = NaN(1, N);            % TR_can;                                               %[mm/s] Transpiration from Vegetation
TR_can_all_store = NaN(1, N, nspecies);                      %[mm/s] Transpiration from Vegetation for each species
Rnrad_can_store = NaN(1, N);            % Rnrad_can;                                         %[W/m^2] Net Radiation Canopy


% Soil Fluxes
Fc_soil_store = NaN(1, N, nspecies);            % VARIABLES.SOIL.Fc_soil;                                             %[umol CO2/ m^2 ground / s] Soil Respiration
H_soil_store = NaN(1, N, nspecies);            % VARIABLES.SOIL.H_soil;                                               %[W/m2] Sensible Heat from Soil Surface
LE_soil_store = NaN(1, N, nspecies);            % VARIABLES.SOIL.LE_soil;                                             %[W/m2] Latent Heat from Soil Surface
G_store = NaN(1, N, nspecies);            % VARIABLES.SOIL.G;                                                         %[W/m2] Ground Heat Flux
Tsurf_store = NaN(1, N, nspecies);                                                 %[C] Soil Surface Temperature
Rnrad_soil_store = NaN(1, N);            % VARIABLES.SOIL.Rnrad_soil;                                       %[W/m2] Net Radiation
RH_soil_store = NaN(1, N);            % RH_soil;                                             %[] Relative Humidity of top soil layer
%E_soil_store(tt) = E_soil;                                                %[]

% Energy Balance Error
Remain_soil_store = NaN(1, N);                                        %[W/m2] Energy Balance Error Soil Surface Energy Balance
Remain_can_store = NaN(1, N);                                          %[W/m2] Energy Balance Error Canopy Energy Balance
Remain_eco_store = NaN(1, N);                                          %[W/m2] Energy Balance Error Ecosystem Energy Balance


% Flux Profiles
An_sun_prof = NaN(nl_can, N);            %(sum(An_sun.*fLAIz,2)).*LAI_sun;                       %[umol CO2/ m^2 ground / s] Photosynthesis minus leaf respiration sunlit fraction for all layers
LE_sun_prof = NaN(nl_can, N);            %(sum(LE_sun.*fLAIz,2)).*LAI_sun;                       %[W/m^2 ground] Latent Heat from sunlit fraction in vegetation for all layers
H_sun_prof = NaN(nl_can, N);            %(sum(H_sun.*fLAIz,2)).*LAI_sun;                         %[W/m^2 ground] Sensible Heat sunlit fraction in vegetation for all layers
Rnrad_sun_prof = NaN(nl_can, N);            %Rnrad_sun;                                          %[W/m^2 ground] Net Radiation Sunlit Fraction in vegetation for all layers
TR_sun_prof = NaN(nl_can, N);            %(sum(TR_sun.*fLAIz,2)).*LAI_sun;                       %[mm/s] Transpiration sunlit fraction in vegetation for all layers
An_shade_prof = NaN(nl_can, N);            %(sum(An_shade.*fLAIz,2)).*LAI_shade;                 %[umol CO2/ m^2 ground / s] Photosynthesis minus leaf respiration shade fraction for all layers
LE_shade_prof = NaN(nl_can, N);            %(sum(LE_shade.*fLAIz,2)).*LAI_shade;                 %[W/m^2 ground] Latent Heat from shade fraction in vegetation for all layers
H_shade_prof = NaN(nl_can, N);            %(sum(H_shade.*fLAIz,2)).*LAI_shade;                   %[W/m^2 ground] Sensible Heat from shade fraction in vegetation for all layers
Rnrad_shade_prof = NaN(nl_can, N);            %Rnrad_shade;                                      %[W/m^2 ground] Net Radiation Shade Fraction in vegetation for all layers
TR_shade_prof = NaN(nl_can, N);            %(sum(TR_shade.*fLAIz,2)).*LAI_sun;                   %[mm/s] Transpiration shade fraction in vegetation for all layers
% Store prof and species for some variables only (per unit of LAI)
Ph_sun_prof_all = NaN(nl_can, N, nspecies);            %VARIABLES.CANOPY.Ph_sun;                         %[umol CO2/ m^2 leaf / s] Photosynthesis sunlit fraction for all layers, and for all species
An_sun_prof_all = NaN(nl_can, N, nspecies);            %VARIABLES.CANOPY.An_sun;                         %[umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration sunlit fraction for all layers, and for all species
Cs_sun_prof_all = NaN(nl_can, N, nspecies);            %VARIABLES.CANOPY.Cs_sun;                         %[umol / mol] Concentration of CO2 in Leaf Boundary Layer sunlit fraction
Ph_shade_prof_all = NaN(nl_can, N, nspecies);            %VARIABLES.CANOPY.Ph_shade;                     %[umol CO2/ m^2 leaf / s] Photosynthesis shade fraction for all layers, and for all species
An_shade_prof_all = NaN(nl_can, N, nspecies);            %VARIABLES.CANOPY.An_shade;                     %[umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration shade fraction for all layers, and for all specie
Cs_shade_prof_all = NaN(nl_can, N, nspecies);            %VARIABLES.CANOPY.Cs_shade;                     %[umol / mol] Concentration of CO2 in Leaf Boundary Layer shade fraction

% Mean Flux Profiles
An_canopy_prof = NaN(nl_can, N);            %An_sun_prof(:,tt) + An_shade_prof(:,tt);            %[umol CO2/ m^2 ground / s] Total Photosynthesis minus leaf respiration for all layers
LE_canopy_prof = NaN(nl_can, N);            %LE_sun_prof(:,tt) + LE_shade_prof(:,tt);            %[W/m^2/ground] Total Latent Heat from Vegetation in all layers
H_canopy_prof = NaN(nl_can, N);            %H_sun_prof(:,tt) + H_shade_prof(:,tt);               %[W/m^2/ground] Total Sensible Heat from Vegetation in all layers
Rnrad_canopy_prof = NaN(nl_can, N);            %Rnrad_sun_prof(:,tt) + Rnrad_shade_prof(:,tt);   %[W/m^2/ground] Net Radiation in the vegetation for all layers
TR_canopy_prof = NaN(nl_can, N);            %TR_sun_prof(:,tt) + TR_shade_prof(:,tt);            %[mm/s] Transpiration for all layers

% Normalized Flux Profiles (ie. per unit LAI)
%   Sunlit
An_sun_norm_prof = NaN(nl_can, N, nspecies);            %An_sun;                                         %[umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration sunlit fraction for all layers
LE_sun_norm_prof = NaN(nl_can, N, nspecies);            %LE_sun;                                         %[W/ m^2 leaf / s] Latent Heat sunlit fraction for all layers
H_sun_norm_prof = NaN(nl_can, N, nspecies);            %H_sun;                                           %[W/ m^2 leaf / s] Sensible Heat sunlit fraction for all layers
Rnrad_sun_norm_prof = NaN(nl_can, N);            %Rnrad_sun./LAI_sun;                            %[W/ m^2 leaf / s] Net Radiation sunlit fraction for all layers

%   Shaded
An_shade_norm_prof = NaN(nl_can, N, nspecies);            %An_shade;                                     %[umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration shade fraction for all layers
LE_shade_norm_prof = NaN(nl_can, N, nspecies);            %LE_shade;                                     %[W/ m^2 leaf / s] Latent Heat shade fraction for all layers
H_shade_norm_prof = NaN(nl_can, N, nspecies);            %H_shade;                                       %[W/ m^2 leaf / s] Sensible Heat shade fraction for all layers
Rnrad_shade_norm_prof = NaN(nl_can, N);            %Rnrad_shade./LAI_shade;                      %[W/ m^2 leaf / s] Net Radiation shade fraction for all layers

%   Canopy
An_canopy_norm_prof = NaN(nl_can, N);            %(sum(An_sun.*fLAIz,2)).*fsun + (sum(An_shade.*fLAIz,2)).*fshade; % [umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration  for all layers
LE_canopy_norm_prof = NaN(nl_can, N);            %(sum(LE_sun.*fLAIz,2)).*fsun + (sum(LE_shade.*fLAIz,2)).*fshade; % [W/ m^2 leaf / s] Latent Heat for all layers
H_canopy_norm_prof = NaN(nl_can, N);            %(sum(H_sun.*fLAIz,2)).*fsun + (sum(H_shade.*fLAIz,2)).*fshade;    % [W/ m^2 leaf / s] Sensible Heat for all layers
Rnrad_canopy_norm_prof = NaN(nl_can, N);            %(Rnrad_sun.*fsun) + (Rnrad_shade.*fshade);                    % [W/ m^2 leaf / s] Net Radiation for all layers

% CANOPY STATES:
% Leaf States
Tl_sun_prof = NaN(nl_can, N, nspecies);            %Tl_sun;                                              %[C] Leaf Temperature Sunlit Fraction all layers and species
Tl_shade_prof = NaN(nl_can, N, nspecies);            %Tl_shade;                                          %[C] Leaf Temperature Shade Fraction all Layers and species
Tl_sun_Ta_Diff = NaN(nl_can, N, nspecies);            %Tl_sun - repmat(TAz,1,nspecies);                  %[C] Temperature Difference Atmosphere and Leaf Sunlint for all layers and species
Tl_shade_Ta_Diff = NaN(nl_can, N, nspecies);            %Tl_shade - repmat(TAz,1,nspecies);              %[C] Temperature Difference Atmosphere and Leaf Shade for all layers and species

psil_sun_prof = NaN(nl_can, N, nspecies);            %psil_sun;                                          %[MPa] Leaf Water Potential sunlit fraction, for all layers and species
psil_shade_prof = NaN(nl_can, N, nspecies);            %psil_shade;                                      %[MPa] Leaf Water Potential Shade fraction, for all layers and species

fsvg_sun_prof = NaN(nl_can, N, nspecies);            %fsvg_sun;                                          %[] g Tuzet factor to Ball Berry Model for sunlit fraction, for all layers and species
fsvm_sun_prof = NaN(nl_can, N, nspecies);            %fsvm_sun;                                          %[] m Tuzet factor to Ball Berry Model for sunlit fraction, for all layers and species
fsvg_shade_prof = NaN(nl_can, N, nspecies);            %fsvg_shade;                                      %[] g Tuzet factor to Ball Berry Model for shade fraction, for all layers and species
fsvm_shade_prof = NaN(nl_can, N, nspecies);            %fsvm_shade;                                      %[] m Tuzet factor to Ball Berry Model for shade fraction, for all layers and species

gsv_sun_prof = NaN(nl_can, N, nspecies);            %gsv_sun;                                            %[mol/m^2/s] Stomatal Conductance for sunlit fraction for all layers and species
gsv_shade_prof = NaN(nl_can, N, nspecies);            %gsv_shade;                                        %[mol/m^2/s] Stomatal Conductance for shade fraction for all layers and species

Ci_sun_prof = NaN(nl_can, N, nspecies);            %Ci_sun;                                              %[umol/mol] Internal leaf concentration of CO2 for sunlit fraction, for all layers and species
Ci_shade_prof = NaN(nl_can, N, nspecies);            %Ci_shade;                                          %[umol/mol] Internal leaf concentration of CO2 for shade fraction, for all layers and species

gbv_sun_prof = NaN(nl_can, N, nspecies);            %gbv_sun;                                            %[mol/m^2/s] Boundary Layer conductance to vapor in sunlit, for all layers and species
gbh_sun_prof = NaN(nl_can, N, nspecies);            %gbh_sun;                                            %[mol/m^2/s] Boundary Layer conductance to heat in sunlit, for all layers and species
gbv_shade_prof = NaN(nl_can, N, nspecies);            %gbv_shade;                                        %[mol/m^2/s] Boundary Layer conductance to vapor in shade, for all layers and species
gbh_shade_prof = NaN(nl_can, N, nspecies);            %gbh_shade;                                        %[mol/m^2/s] Boundary Layer conductance to heat in shade, for all layers and species

LAI_sun_prof = NaN(nl_can, N);            %LAI_sun;                                              %[m^2 leaf / m^2 ground] Leaf Area Index sunlit, for all layers
LAI_shade_prof = NaN(nl_can, N);            %LAI_shade;                                          %[m^2 leaf / m^2 ground] Leaf Area Index shade, for all layers

fsun_prof = NaN(nl_can, N);            %fsun;                                                    %[] Fraction of LAI sunlit
fshade_prof = NaN(nl_can, N);            %fshade;                                                %[] Fraction of LAI shade

LSsunCON_store = NaN(1, N, nspecies);                                           % Convergence of leaf solution, for sunlit, for all species
LSshaCON_store = NaN(1, N, nspecies);                                          % Convergence of leaf solution, for shade, for all specie

% leaf convergence
cntleafsolution_store = NaN(1, N, nspecies);            % VARIABLES.CANOPY.cntleafsolution;          % Leaf convergence

% Photosynthetic Biochemistry
Ph_limit_sun_store = NaN(nl_can, N, nspecies);            %Ph_limit_sun;                                 %[] Type of Photosynthesis 1. Rubisco-limited, 2light-limited, 3. Sucrose Limited
Jc_C3_sun_store = NaN(nl_can, N, nspecies);            %Jc_C3_sun;                                       %[umol/m^2 leaf area/s] Rubisco-Limited Rate for sunlit, for all layers, and species
Jj_C3_sun_store = NaN(nl_can, N, nspecies);            %Jj_C3_sun;                                       %[umol/m^2 leaf area/s] Light-Limited Rate for sunlit, for all layers, and species
Js_C3_sun_store = NaN(nl_can, N, nspecies);            %Js_C3_sun;                                       %[umol/m^2 leaf area/s] Sucrose-Limited Rate for sunlit, for all layers, and species
Jc_C4_sun_store = NaN(nl_can, N, nspecies);            %Jc_C4_sun;                                       %[umol/m^2 leaf area/s] Rubisco-Limited Rate for sunlit, for all layers, and species
Jj_C4_sun_store = NaN(nl_can, N, nspecies);            %Jj_C4_sun;                                       %[umol/m^2 leaf area/s] Ligth-Limited Rate for sunlit, for all layers, and species
Js_C4_sun_store = NaN(nl_can, N, nspecies);            %Js_C4_sun;                                       %[umol/m^2 leaf area/s] Sucrose-Limited Rate for sunlit, for all layers, and species

Ph_limit_shade_store = NaN(nl_can, N, nspecies);            %Ph_limit_shade;                             %[] Type of Photosynthesis 1. Rubisco-limited, 2light-limited, 3. Sucrose Limited
Jc_C3_shade_store = NaN(nl_can, N, nspecies);            %Jc_C3_shade;                                   %[umol/m^2 leaf area/s] Rubisco-Limited Rate for shade, for all layers, and species
Jj_C3_shade_store = NaN(nl_can, N, nspecies);            %Jj_C3_shade;                                   %[umol/m^2 leaf area/s] Ligth-Limited Rate for shade, for all layers, and species
Js_C3_shade_store = NaN(nl_can, N, nspecies);            %Js_C3_shade;                                   %[umol/m^2 leaf area/s] Sucrose-Limited Rate for shade, for all layers, and species
Jc_C4_shade_store = NaN(nl_can, N, nspecies);            %Jc_C4_shade;                                   %[umol/m^2 leaf area/s] Rubisco-Limited Rate for shade, for all layers, and species
Jj_C4_shade_store = NaN(nl_can, N, nspecies);            %Jj_C4_shade;                                   %[umol/m^2 leaf area/s] Ligth-Limited Rate for shade, for all layers, and species
Js_C4_shade_store = NaN(nl_can, N, nspecies);            %Js_C4_shade;                                   %[umol/m^2 leaf area/s] Sucrose-Limited Rate for shade, for all layers, and species

dryfrac_prof = NaN(nl_can, N);            %dryfrac;                                              %[] Fraction of LAI that is dry
wetfrac_prof = NaN(nl_can, N);            %wetfrac;                                              %[] Fraction of LAI that is wet

ppt_ground_store = NaN(1, N, nspecies);            %ppt_ground;                                       %[mm] Precipitation Reaching the ground

% Mean Canopy Profiles
Ci_canopy_prof = NaN(nl_can, N);            %(nansum(Ci_sun.*fLAIz,2)).*fsun + (nansum(Ci_shade.*fLAIz,2)).*fshade;       %[umol/mol] Internal leaf concentration of CO2, for all layers
Tl_canopy_prof = NaN(nl_can, N);            %(nansum(Tl_sun.*fLAIz,2)).*fsun + (nansum(Tl_shade.*fLAIz,2)).*fshade;       %[C] Leaf Temperature, for all layers
gsv_canopy_prof = NaN(nl_can, N);            %(nansum(gsv_sun.*fLAIz,2)).*fsun + (nansum(gsv_shade.*fLAIz,2)).*fshade;    %[mol/m^2/s] Stomatal Conductance for all layers
psil_canopy_prof = NaN(nl_can, N);            %(nansum(psil_sun.*fLAIz,2)).*fsun + (nansum(psil_shade.*fLAIz,2)).*fshade; %[MPa] Leaf Water Potential for all layers
fsvg_canopy_prof = NaN(nl_can, N);            %(nansum(fsvg_sun.*fLAIz,2)).*fsun + (nansum(fsvg_shade.*fLAIz,2)).*fshade; %[] g Tuzet factor to Ball Berry Model for all layers and species
fsvm_canopy_prof = NaN(nl_can, N);            %(nansum(fsvm_sun.*fLAIz,2)).*fsun + (nansum(fsvm_shade.*fLAIz,2)).*fshade; %[] m Tuzet factor to Ball Berry Model for all layers and species
gbv_canopy_prof = NaN(nl_can, N);            %(nansum(gbv_sun.*fLAIz,2)).*fsun + (nansum(gbv_shade.*fLAIz,2)).*fshade;    %[mol/m^2/s] Boundary Layer conductance to vapor for all layers
gbh_canopy_prof = NaN(nl_can, N);            %(nansum(gbh_sun.*fLAIz,2)).*fsun + (nansum(gbh_shade.*fLAIz,2)).*fshade;    %[mol/m^2/s] Boundary Layer conductance to heat for all layers

% Mean Canopy States
Ci_mean = NaN(1, N);            % sum(Ci_canopy_prof(vinds,tt).*LADnorm(vinds))/sum(LADnorm);                         %[umol/mol] Internal leaf concentration of CO2, mean canopy
Tl_mean = NaN(1, N);            % sum(Tl_canopy_prof(vinds,tt).*LADnorm(vinds))/sum(LADnorm);                         %[C] Leaf Temperature of CO2, mean canopy
gsv_mean = NaN(1, N);            % sum(gsv_canopy_prof(vinds,tt).*LADnorm(vinds))/sum(LADnorm);                       %[mol/m^2/s] Stomatal Conductance, mean canopy
psil_mean = NaN(1, N);            % sum(psil_canopy_prof(vinds,tt).*LADnorm(vinds))/sum(LADnorm);                     %[MPa] Leaf Water Potential for all layers
fsvg_mean = NaN(1, N);            % sum(fsvg_canopy_prof(vinds,tt).*LADnorm(vinds))/sum(LADnorm);                     %[] g Tuzet factor to Ball Berry Model for all layers and specie
fsvm_mean = NaN(1, N);            % sum(fsvm_canopy_prof(vinds,tt).*LADnorm(vinds))/sum(LADnorm);                     %[] m Tuzet factor to Ball Berry Model for all layers and specie


% Canopy Microenvironment
CAz_prof = NaN(nl_can, N, nspecies);            %VARIABLES.CANOPY.CAz;                                                      %[umol/mol] Atmospheric Concentration of CO2, for all layers
TAz_prof = NaN(nl_can, N, nspecies);            %VARIABLES.CANOPY.TAz;                                                      %[C] Air Temperature, for all layers
EAz_prof = NaN(nl_can, N, nspecies);            %VARIABLES.CANOPY.EAz;                                                      %[kPa] Vapor Pressure Air, for all layers
Uz_prof = NaN(nl_can, N);            %VARIABLES.CANOPY.Uz;                                                        %[m/s] Wind Velocity

% Canopy H2O Storage
Sh2o_canopy_prof = NaN(nl_can, N);            %Sh2o_prof;                                        %[mm] Intercepted water in canopy, for all layers
Sh2o_canopy = NaN(1, N);            % Sh2o_can;                                              %[mm] Total Intercepted water in canopy.
% VARIABLES.CANOPY.Sh2o_can_prev = Sh2o_can;                                 %[mm] Total Intercepted water in canopy in previous time step

% Condensation
Ch2o_canopy_prof = NaN(nl_can, N);            %Ch2o_prof;                                        %[mm] Condensation water in the canopy, for all layers
Ch2o_canopy = NaN(1, N);            % Ch2o_can;                                              %[mm] Total Condensation water in the canopy.

% Evaporation
Evap_canopy_prof = NaN(nl_can, N);            %Evap_prof;                                        %[mm] Evaporation water from canopy
Evap_canopy = NaN(1, N);            % Evap_can;                                              %[mm] Total Evaporation water from canopy.


% SOIL VARIABLES
volliq_store = NaN(nl_soil, N, nspecies);            %volliq;                                               %[] volumetric water content
krad_store = NaN(nl_soil, N, nspecies);            %krad;                                                 %[/s] radial conductivity of the root
kax_store = NaN(nl_soil, N, nspecies);            %kax;                                                   %[mm/s] axial conductivity of the root
hk_store = NaN(nl_soil, N);            %hk;                                                       %[mm/s] Soil Hydraulic Conductivity at Node Layers
smp_store = NaN(nl_soil, N, nspecies);            %smp;                                                     %[mm of H2O] Soil Matric Potential
rpp_store = NaN(nl_soil, N, nspecies);            %rpp;                                                   %[mm of H2O] Root Water Potential
mberrormm_store = NaN(1, N);                                           %[mm/dtime] Mass balance error.


smpMPA_store = NaN(nl_soil, N, nspecies);            %smp*CONSTANTS.mmH2OtoMPa;                             %[MPa] Soil Matric Potential for all layers
rppMPA_store = NaN(nl_soil, N, nspecies);            %rpp*CONSTANTS.mmH2OtoMPa;                           %[MPa] Root Water Potential for all layers and species

% smp_weight_store(tt) = smp_wgt*CONSTANTS.mmH2OtoMP                     %[MPA] Weighted mean smp over root uptake profile
rpp_weight_store = NaN(1, N, nspecies);            % rpp_wgt*CONSTANTS.mmH2OtoMPa;                %[Mpa] Weighted mean Root weight to use in Canopy Model
%    volliq_weight_store(tt) = sum(volliq .* VERTSTRUC.rootfr);

qlayer_store = NaN(nl_soil+1, N);            %VARIABLES.SOIL.qlayer;                                %[mm/s] Fluxes of water between layers in the soil
bot_drainage_store = NaN(1, N);            % VARIABLES.SOIL.qlayer(end);                                    %[mm/s] Flux of water at bottom of soil column
hor_drainage_store = NaN(1, N);            % VARIABLES.SOIL.hor_drainage;                    %[mm/s] Horizontal Drainage. Only happens when superaturation occurs.
hor_drainage_lay_store = NaN(nl_soil, N);            %VARIABLES.SOIL.hor_drainage_lay;            %[mm/s] Horizontal Drainage for all layers. Only happens when superaturation occurs.


if (SWITCHES.soilheat_on)
    Ts_store = NaN(nl_soil, N, nspecies);            %VARIABLES.SOIL.Ts;                                %[C] Temperature in the soil. For all layers
    Gnew_store = NaN(nl_soil, N, nspecies);            %VARIABLES.SOIL.Gnew;                            %[C] Soil Heat Flux. For all layers in the soil.
    cpv_store = NaN(nl_soil, N, nspecies);            %VARIABLES.SOIL.cpv;                              %[J / m^3 / K] Soil Volumetric Heat Capacity
    TKsoil_store = NaN(nl_soil, N, nspecies);            %VARIABLES.SOIL.TKsoil;                        %[W / m / K] Thermal Conductivity Soil
    TKsol_store = NaN(nl_soil, N, nspecies);            %VARIABLES.SOIL.TKsol;                          %[W / m / K] Thermal Conductivity of Mineral Solids
end

if (SWITCHES.ns)
    wuptake_store = NaN(nl_soil, N);            %layeruptake;                                 %[mm/s] Fluxes of plant water uptake at each layer. Shows the total from all the species
    wuptake_all_store = NaN(nl_soil, N, nspecies);            %layeruptake_all;                        %[mm/s] Fluxes of plant water uptake at each layer, and for all species.
end


runoff_store = NaN(1, N, nspecies);            % VARIABLES.SOIL.runoff;                                %[mm] Flux of water lost as runoff
qback_store = NaN(1, N, nspecies);            % VARIABLES.SOIL.qback;                                  %[mm/s] Flux of water returned to litter or surface by sypersaturation of first layer
qadflux_store = NaN(1, N, nspecies);            % VARIABLES.SOIL.qadflux;                              %[mm/s]  Final Flux that Infiltrates in the soil
qback_lit_store = NaN(1, N, nspecies);            % VARIABLES.SOIL.qback_lit;                          %[mm/s] Real flux that is going back to snow-litter pack

Gs_store = NaN(1, N, nspecies);            % VARIABLES.SOIL.Gs;                                        %Ground Heat Flux into the Soil from the Soil- Snow-Litter Pack  Boundary [W/m2]
Gsl_store = NaN(1, N, nspecies);          
LE_soli_store = NaN(1, N, nspecies);      
LE_liat_store = NaN(1, N, nspecies);      
psili_MPa_store = NaN(1, N, nspecies);    
psili_store = NaN(1, N, nspecies);        
Tli_store = NaN(1, N, nspecies);          
Tsl_store = NaN(1, N, nspecies);          
% rhosn
rhosn_store = NaN(1, N, nspecies);        
% w
wliqsl_store = NaN(1, N, nspecies);       
wicesl_store = NaN(1, N, nspecies);       
deltawice_store = NaN(1, N, nspecies);    
dH_store = NaN(1, N, nspecies);           
dS_store = NaN(1, N, nspecies);           

wsn_store = NaN(1, N, nspecies);          
% z
zicesl_store = NaN(1, N, nspecies);        
zliqsl_store = NaN(1, N, nspecies);        

zsn_store = NaN(1, N, nspecies);           
% vols
volliqli_store = NaN(1, N, nspecies);          
voliceli_store = NaN(1, N, nspecies);          
volliqsn_store = NaN(1, N, nspecies);         
volicesn_store = NaN(1, N, nspecies);         
voltotsn_store = NaN(1, N, nspecies);         
voltotli_store = NaN(1, N, nspecies);         
%Other variables
CRm_store = NaN(1, N, nspecies);          
CRo_store = NaN(1, N, nspecies);          
drhosn_store = NaN(1, N, nspecies);       
qinflL_store = NaN(1, N, nspecies);       
net_qinflL_store = NaN(1, N, nspecies);   
qinfl_store = NaN(1, N, nspecies);        
drainlitter_store = NaN(1, N, nspecies);  
pptrate_ground_store = NaN(1, N, nspecies);
Esoil_store = NaN(1, N, nspecies);         
Esl_store = NaN(1, N, nspecies);           
snow_tcount_store = NaN(1, N, nspecies);        


%% Rohit start
% CROP GROWTH

LeafCMass_can_store = NaN(1, N, nspecies);       
StemCMass_can_store = NaN(1, N, nspecies);       
RootCMass_can_store = NaN(1, N, nspecies);       
GrainCMass_can_store = NaN(1, N, nspecies);      

LAIsim_can_store = NaN(1, N, nspecies);           

TotalCMass_can_store = NaN(1, N, nspecies);       
AboveCMass_can_store = NaN(1, N, nspecies);       


LeafBiomass_can_store = NaN(1, N, nspecies);     
StemBiomass_can_store = NaN(1, N, nspecies);     
RootBiomass_can_store = NaN(1, N, nspecies);     
GrainBiomass_can_store = NaN(1, N, nspecies);    

TotalBioMass_can_store = NaN(1, N, nspecies);    
AboveBioMass_can_store = NaN(1, N, nspecies);    

GrainYield_store = NaN(1, N, nspecies);          

DeathRate_store = NaN(1, N, nspecies);      
NDemand_store = NaN(1, N, nspecies);        

Nreq_store = NaN(1, N, nspecies);       
Fert_req_store = NaN(1, N, nspecies);            
Fert_doy_store = NaN(1, N, nspecies);            

K_ec_store = NaN(1, N, nspecies);
K_rub_leaf_store = NaN(1, N, nspecies);
K_lhc_leaf_store = NaN(1, N, nspecies);
fn_store = NaN(1, N, nspecies);

% Irrigation

IW_store = NaN(1, N, nspecies);          
% 
% VARIABLES.CANOPY.psil_mean = sum(psil_canopy_prof(vinds,tt).*LADnorm(vinds))/sum(LADnorm);
% VARIABLES.CANOPY.gsv_mean = sum(gsv_canopy_prof(vinds,tt).*LADnorm(vinds))/sum(LADnorm);
% %
% Tl_canopy_prof_up = NaN(nl_can, N);            %(nansum(Tl_sun_up.*fLAIz,2)).*fsun + (nansum(Tl_shade_up.*fLAIz,2)).*fshade;
%
% Tl_canopy_prof_lo = NaN(nl_can, N);            %(nansum(Tl_sun_lo.*fLAIz,2)).*fsun + (nansum(Tl_shade_lo.*fLAIz,2)).*fshade;
%
%
Tl_mean_up_store = NaN(1, N, nspecies);       
Tl_mean_lo_store = NaN(1, N, nspecies);      

% VARIABLES.CANOPY.Tl_mean_up = sum(Tl_canopy_prof_up(vinds,tt).*LADnorm(vinds))/sum(LADnorm);
%
% VARIABLES.CANOPY.Tl_mean_lo = sum(Tl_canopy_prof_lo(vinds,tt).*LADnorm(vinds))/sum(LADnorm);

CWSI_store = NaN(1, N, nspecies);           

% Nitrogen Canopy

% Vcmax_store = NaN(1, N);            
Vcmax_vert_store = NaN(nl_can, N, nspecies);            
Chl_la_store = NaN(1, N, nspecies);            
N_can_store = NaN(1, N, nspecies);            
N_leaf_la_store = NaN(nl_can, N, nspecies);   
N_rub_la_store = NaN(nl_can, N, nspecies);    

N_uptake_can_store = NaN(1, N, nspecies);     

Nl_ncan_ga_store = NaN(1, N, nspecies);       
Nr_ncan_ga_store = NaN(1, N, nspecies);       
Ns_ncan_ga_store = NaN(1, N, nspecies);       
Ng_ncan_ga_store = NaN(1, N, nspecies);       


% Mass Balance


MBerror_soil= NaN(1, N, nspecies);
MBerror_littersoil= NaN(1, N, nspecies);
MBerror_mbcan= NaN(1, N, nspecies);
MBerror_mbcanlittersoil= NaN(1, N, nspecies);



MB_store.PPTf = nan(1,N);                   % Total Rainfall  [mm/s]
MB_store.PPTgr = nan(1,N);                  % Total Rainfall Reaching Ground  [mm/s]
MB_store.INFsoil = nan(1,N);                % Infiltration into Soil, Drainage From Litter [mm/s]
MB_store.EVli = nan(1,N);                   % Evaporation from litter [mm/s]
MB_store.EVsoil = nan(1,N);                 % Evaporation from soil [mm/s]
MB_store.EVcanf = nan(1,N);                 % Evaporation from canopy [mm/s]
MB_store.drainage = nan(1,N);               % Bottom Drainage from soil [mm/s]
MB_store.TR = nan(1,N,nspecies);            % Transpiration all species [mm/s]
MB_store.runoff = nan(1,N);                 % Transpiration all species [mm/s]
MB_store.COcanf = nan(1,N);                 % Total canopy condensation [mm/s]

MB_store.dws_f = nan(1,N);                  % Rate of change in soil water content [mm/s]
MB_store.dwc_f = nan(1,N);                  % Rate of change in canopy water content [mm/s]
MB_store.dwl_f = nan(1,N);                  % Rate of change in litter water content [mm/s]

MB_store.mbsoil = nan(1,N);                 % rate of change in soil mass balance content[mm/day]
MB_store.mblittersoil = nan(1,N);           % rate of change in soil - litter mass balance content [mm/day]
MB_store.mbcan = nan(1,N);                  % rate of change in canopy mass balance content [mm/day]
MB_store.mbcanlittersoil = nan(1,N);        % rate of change in soil - litter - canopy mass balance content [mm/day]



if (SWITCHES.entropy_on)
    % CANOPY
    % in
    SScandir_in_all_lay_store = nan(nl_can,N,nspecies);
    SScandir_in_all_store = nan(1,N,nspecies);
    SScandir_in_tot_store = nan(1,N);
    
    SScandif_in_all_lay_store = nan(nl_can,N,nspecies) ;
    SScandif_in_all_store = nan(1,N,nspecies);
    SScandif_in_tot_store = nan(1,N);
    
    SScanLW_in_all_lay_store = nan(nl_can,N,nspecies);
    SScanLW_in_all_store = nan(1,N,nspecies);
    SScanLW_in_tot_store = nan(1,N);
    
    SScan_in_all_lay_store = nan(nl_can,N,nspecies);
    SScan_in_all_store = nan(1,N,nspecies);
    SScan_in_store = nan(1,N);
    
    % out
    SScandir_out_all_lay_store = nan(nl_can,N,nspecies);
    SScandir_out_all_store = nan(1,N,nspecies);
    SScandir_out_tot_store = nan(1,N);
    
    SScandif_out_all_lay_store = nan(nl_can,N,nspecies);
    SScandif_out_all_store = nan(1,N,nspecies);
    SScandif_out_tot_store = nan(1,N);
    
    SScanLW_out_all_lay_store = nan(nl_can,N,nspecies);
    SScanLW_out_all_store = nan(1,N,nspecies);
    SScanLW_out_tot_store = nan(1,N);
    
    SScanLE_out_all_lay_store = nan(nl_can,N,nspecies);
    SScanLE_out_all_store = nan(1,N,nspecies);
    SScanLE_out_tot_store = nan(1,N);
    
    SScanH_out_all_lay_store = nan(nl_can,N,nspecies);
    SScanH_out_all_store = nan(1,N,nspecies);
    SScanH_out_tot_store = nan(1,N);
    
    % Total OUT
    SScan_out_all_lay_store = nan(nl_can,N,nspecies);
    SScan_out_all_store = nan(1,N,nspecies);
    SScan_out_store = nan(1,N);
    
    % TOTAL CANOPY
    SScan_all_lay_store = nan(nl_can,N,nspecies);
    SScan_all_store = nan(1,N,nspecies);
    SScan_tot_store = nan(1,N);
    
    % PHOTOSYNTHESIS
    Epho_dif_store = nan(1,N);
    Epho_dir_store = nan(1,N);
    Epho_tot_store = nan(1,N);
    SSphodir_in_store = nan(1,N);
    SSphodif_in_store = nan(1,N);
    
    
    % SOIL
    % in
    SSsoildir_in_tot_store = nan(1,N);
    SSsoildif_in_tot_store = nan(1,N);
    SSsoilLW_in_tot_store = nan(1,N);
    SSsoilG_in_tot_store = nan(1,N);
    SSsoil_in_store = nan(1,N);
    
    % out
    SSsoildir_out_tot_store = nan(1,N);
    SSsoildif_out_tot_store = nan(1,N);
    SSsoilLW_out_tot_store = nan(1,N);
    SSsoilLE_out_tot_store = nan(1,N);
    SSsoilH_out_tot_store = nan(1,N);
    SSsoil_out_store = nan(1,N);
    
    % G, dS, dH SOIL
    SSsoilG_store = nan(1,N);
    SSsoildH_store = nan(1,N);
    SSsoildS_store = nan(1,N);
    
    % dH CANOPY
    SSdHcan_store = nan(1,N);
    
    % TOTAL SOIL
    SSsoil_tot_store = nan(1,N);
    
    
    % TOTALS BY TYPE OF ENERGY
    SSeco_totSWdir_store = nan(1,N);
    SSeco_totSWdif_store = nan(1,N);
    SSeco_totLW_store = nan(1,N);
    SSeco_totH_store = nan(1,N);
    SSeco_totLE_store = nan(1,N);
    SSeco_tot_store_test = nan(1,N);
    
    
    % TOTAL ECOSYSTEM
    SSeco_tot_store = nan(1,N);
    
    % NET CANOPY
    
    SSdir_net_in_store = nan(1,N);
    SSdir_net_out_store = nan(1,N);
    SSdif_net_in_store = nan(1,N);
    SSdif_net_out_store = nan(1,N);
    
    SSLW_net_in_store = nan(1,N);
    SSLW_net_inX_store = nan(1,N);
    SSLW_net_out_store = nan(1,N);
    SSLW_net_outX_store = nan(1,N);
    
    SSSW_net_inX_store = nan(1,N);
    SSSW_net_outX_store = nan(1,N);
    
    
    SSLE_net_store = nan(1,N);
    SSH_net_store = nan(1,N);
    SSG_net_store = nan(1,N);
    SSdHcan_net_store = nan(1,N);
    
    SSeco_net_in_store = nan(1,N);
    SSeco_net_out_store = nan(1,N);
    SSeco_net_tot_store = nan(1,N);
    SSeco_net_ratio_store = nan(1,N);
    
    SSTl_net_store = nan(1,N);
    SSTl_net2_store = nan(1,N);
    
    % ENERGY AND TEFF
    
    SSeco_tot_store = nan(1,N);
    SWdir_in_store = nan(1,N);
    SWdir_out_store = nan(1,N);
    SWdif_in_store = nan(1,N);
    SWdif_out_store = nan(1,N);
    LWin_net_store = nan(1,N);
    LWout_net_store = nan(1,N);
    LWemi_net_store = nan(1,N);
    LE_net_store = nan(1,N);
    H_net_store = nan(1,N);
    G_net_store = nan(1,N);
    dHcan_net_store = nan(1,N);
    SSTeffent_store = nan(1,N);
    SSXeffent_store = nan(1,N);
    
    SH2O_store = nan(nl_soil-1,N);
    EH2O_store = nan(nl_soil-1,N);
    
end


