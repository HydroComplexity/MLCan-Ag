%=========================================================================
% This script stores some of the variables in the code
%
% Written by Darren Drewry and Juan Camilo Quijano, UIUC, 2013


tpq=tt-NoDay_prev*(24*60/CONSTANTS.timestep);

% Store simulation results

LAIz = VERTSTRUC.LAIz;                                                      % Leaf area index distributed at each layer in z

zmm=VERTSTRUC.znsmm;                                                        % [mm] Vector with depth at nodes of each layer
dzmm=VERTSTRUC.dzsmm;                                                       % [mm] Vector with Layer thicknesses
zimm=VERTSTRUC.zhsmm;                                                       % [mm] Vector with depth of interfaces between the layers

% CANOPY VARIABLES
% volliq_store(:,tpq,:) = volliq;                                               %[] volumetric water content  

% Radiation Absorption
PARabs_sun_prof(:,tpq) = PARabs_sun;                                        %[umol / m^2 ground / s] absorbed PAR by sunlit fraction in each layer  
PARabs_shade_prof(:,tpq) = PARabs_shade;                                    %[umol / m^2 ground / s] absorbed PAR by shade fraction in each layer%#ok<AGROW>
PARabs_canopy_prof(:,tpq) = PARabs_sun + PARabs_shade;                      %[umol / m^2 ground / s] total absorbed PAR in each layer
                
PARabs_sun_norm_prof(:,tpq) = PARabs_sun./LAIz;                             %[umol / m^2 leaf / s ] absorbed PAR by sunlit fraction normalized by LAI  in each layer
PARabs_shade_norm_prof(:,tpq) = PARabs_shade./LAIz;                         %[umol / m^2 leaf / s ] absorbed PAR by shade fraction normalized by LAI  in each layer
PARabs_canopy_norm_prof(:,tpq) = (PARabs_sun + PARabs_shade)./LAIz;         %[umol / m^2 leaf / s ] absorbed PAR normalized by LAI  in each layer
                
NIRabs_sun_prof(:,tpq) = NIRabs_sun;                                        %[W / m^2 ground] absorbed NIR by sunlit fraction in each layer
NIRabs_shade_prof(:,tpq) = NIRabs_shade;                                    %[W / m^2 ground] absorbed NIR by shade fraction in each layer
NIRabs_canopy_prof(:,tpq) = NIRabs_sun + NIRabs_shade;                      %[W / m^2 ground] total absorbed NIR in each layer
        
NIRabs_sun_norm_prof(:,tpq) = NIRabs_sun./LAIz;                             %[W / m^2 leaf] absorbed NIR by sunlit fraction normalized by LAI in each layer
NIRabs_shade_norm_prof(:,tpq) = NIRabs_shade./LAIz;                         %[W / m^2 leaf] absorbed NIR by shade fraction normalized by LAI in each layer 
NIRabs_canopy_norm_prof(:,tpq) = (NIRabs_sun + NIRabs_shade)./LAIz;         %[W / m^2 leaf] total absorbed NIR normalized by LAI in each layer
        
LWabs_can_prof(:,tpq) = LWabs_can;                                          %[W / m^2 ground] absorbed LW by canopy in each layer
LWemit_can_prof(:,tpq) = LWemit_can;                                        %[W / m^2 ground] emitted LW by canopy in each layer 

LWemit_can_store(1,tpq) = sum(LWabs_can);                                   %[W / m^2 ground] total absorbed LW by canopy        
LWemit_soil_store(1,tpq) = LWemit_soil;                                     %[W / m^2 ground] emitted LW by canopy 
        
SWout_store(1,tpq) = SWout;                                                 %[W / m^2 ground] total reflected SW from canopy
LWout_store(1,tpq) = LWout;                                                 %[W / m^2 ground] total leaving LW from canopy
LWdn_in(1,tpq) = FORCING.LWdn;                                              %[W / m^2 ground] total incoming LW radiation  
fdiff_store(1,tpq) = fdiff;                                                 %[] Fraction of incoming radiation that is diffuse     
refl_soil_store(1,tpq) = VARIABLES.refl_soil;                               %[] Reflected SW radiation by Soil 
        
% FLUXES:       
% Ecosystem Fluxes (Canopy + Soil)    
Fc_eco_store(1,tpq) = Fc_soil - An_can;                                     %[umol CO2/ m^2 ground / s] flux of CO2 from surface, similat to NEE 
LE_eco_store(1,tpq) = LE_can + LE_soil;                                     %[W/m^2] Total latent heat released from the surface  
H_eco_store(1,tpq) = H_can + H_soil;                                        %[W/m^2] Total sensible heat released from the surface
Rnrad_eco_store(1,tpq) = Rnrad_eco;                                         %[W/m^2] Ecosystem Net Radiation
        
% Canopy Fluxes
Ph_can_store(1,tpq) = Ph_can;                                               %[umol CO2/ m^2 ground / s] Photosynthetic flux from ecosystem    
Ph_can_all_store(1,tpq,:)=VARIABLES.CANOPY.Ph_can_all;                      %[umol CO2/ m^2 ground / s] Photosynthetic flux for each species  
An_can_store(1,tpq) = An_can;                                               %[umol CO2/ m^2 ground / s] Photosynthetic flux minus leaf respiration from ecosystem  
An_can_all_store(1,tpq,:)=VARIABLES.CANOPY.An_can_all;                      %[umol CO2/ m^2 ground / s] Photosynthetic flux minus from each species  
LE_can_store(1,tpq) = LE_can;                                               %[W/m^2] Total Latent Heat from Vegetation 
LE_can_all_store(1,tpq,:)=VARIABLES.CANOPY.LE_can_all;                      %[W/m^2] Total Latent Heat from Vegetation for each species   
H_can_store(1,tpq) = H_can;                                                 %[W/m^2] Total Sensible Heat from Vegetation
H_can_all_store(1,tpq,:)=VARIABLES.CANOPY.H_can_all;                        %[W/m^2] Total Sensible Heat from Vegetation for each species
dH_can_store(1,tpq) = dHcan;                                                %[W/m^2] Total change of heat content in leaves
dH_can_all_store(1,tpq,:)=VARIABLES.CANOPY.dHcan_all;                       %[W/m^2] Total change of heat content in leaves for each species 
TR_can_store(1,tpq) = TR_can;                                               %[mm/s] Transpiration from Vegetation  
TR_can_all_store(1,tpq,:)=VARIABLES.CANOPY.TR_can_all;                      %[mm/s] Transpiration from Vegetation for each species  
Rnrad_can_store(1,tpq) = Rnrad_can;                                         %[W/m^2] Net Radiation Canopy  
                
        
   % Soil Fluxes
Fc_soil_store(1,tpq,:) = VARIABLES.SOIL.Fc_soil;                                             %[umol CO2/ m^2 ground / s] Soil Respiration
H_soil_store(1,tpq,:) = VARIABLES.SOIL.H_soil;                                               %[W/m2] Sensible Heat from Soil Surface   
LE_soil_store(1,tpq,:) = VARIABLES.SOIL.LE_soil;                                             %[W/m2] Latent Heat from Soil Surface 
G_store(1,tpq,:) = VARIABLES.SOIL.G;                                                         %[W/m2] Ground Heat Flux
Tsurf_store(1,tpq,:)=VARIABLES.SOIL.Tsurf;                                                   %[C] Soil Surface Temperature  
Rnrad_soil_store(1,tpq) = VARIABLES.CANOPY.Rnrad_soil;                                       %[W/m2] Net Radiation
RH_soil_store(1,tpq) = RH_soil;                                             %[] Relative Humidity of top soil layer 
%E_soil_store(tpq) = E_soil;                                                %[]                                           
   
% Energy Balance Error 
Remain_soil_store(1,tpq)=remainsoil;                                        %[W/m2] Energy Balance Error Soil Surface Energy Balance
Remain_can_store(1,tpq)=remaincan;                                          %[W/m2] Energy Balance Error Canopy Energy Balance
Remain_eco_store(1,tpq)=remaineco;                                          %[W/m2] Energy Balance Error Ecosystem Energy Balance
                        
        
    % Flux Profiles
An_sun_prof(:,tpq) = (sum(An_sun.*fLAIz,2)).*LAI_sun;                       %[umol CO2/ m^2 ground / s] Photosynthesis minus leaf respiration sunlit fraction for all layers
LE_sun_prof(:,tpq) = (sum(LE_sun.*fLAIz,2)).*LAI_sun;                       %[W/m^2 ground] Latent Heat from sunlit fraction in vegetation for all layers
H_sun_prof(:,tpq) = (sum(H_sun.*fLAIz,2)).*LAI_sun;                         %[W/m^2 ground] Sensible Heat sunlit fraction in vegetation for all layers
Rnrad_sun_prof(:,tpq) = Rnrad_sun;                                          %[W/m^2 ground] Net Radiation Sunlit Fraction in vegetation for all layers
TR_sun_prof(:,tpq) = (sum(TR_sun.*fLAIz,2)).*LAI_sun;                       %[mm/s] Transpiration sunlit fraction in vegetation for all layers 
An_shade_prof(:,tpq) = (sum(An_shade.*fLAIz,2)).*LAI_shade;                 %[umol CO2/ m^2 ground / s] Photosynthesis minus leaf respiration shade fraction for all layers
LE_shade_prof(:,tpq) = (sum(LE_shade.*fLAIz,2)).*LAI_shade;                 %[W/m^2 ground] Latent Heat from shade fraction in vegetation for all layers
H_shade_prof(:,tpq) = (sum(H_shade.*fLAIz,2)).*LAI_shade;                   %[W/m^2 ground] Sensible Heat from shade fraction in vegetation for all layers
Rnrad_shade_prof(:,tpq) = Rnrad_shade;                                      %[W/m^2 ground] Net Radiation Shade Fraction in vegetation for all layers        
TR_shade_prof(:,tpq) = (sum(TR_shade.*fLAIz,2)).*LAI_sun;                   %[mm/s] Transpiration shade fraction in vegetation for all layers
        % Store prof and species for some variables only (per unit of LAI)
Ph_sun_prof_all(:,tpq,:) = VARIABLES.CANOPY.Ph_sun;                         %[umol CO2/ m^2 leaf / s] Photosynthesis sunlit fraction for all layers, and for all species
An_sun_prof_all(:,tpq,:) = VARIABLES.CANOPY.An_sun;                         %[umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration sunlit fraction for all layers, and for all species
Cs_sun_prof_all(:,tpq,:) = VARIABLES.CANOPY.Cs_sun;                         %[umol / mol] Concentration of CO2 in Leaf Boundary Layer sunlit fraction
Ph_shade_prof_all(:,tpq,:) = VARIABLES.CANOPY.Ph_shade;                     %[umol CO2/ m^2 leaf / s] Photosynthesis shade fraction for all layers, and for all species 
An_shade_prof_all(:,tpq,:) = VARIABLES.CANOPY.An_shade;                     %[umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration shade fraction for all layers, and for all specie 
Cs_shade_prof_all(:,tpq,:) = VARIABLES.CANOPY.Cs_shade;                     %[umol / mol] Concentration of CO2 in Leaf Boundary Layer shade fraction
               
    % Mean Flux Profiles
An_canopy_prof(:,tpq) = An_sun_prof(:,tpq) + An_shade_prof(:,tpq);            %[umol CO2/ m^2 ground / s] Total Photosynthesis minus leaf respiration for all layers 
LE_canopy_prof(:,tpq) = LE_sun_prof(:,tpq) + LE_shade_prof(:,tpq);            %[W/m^2/ground] Total Latent Heat from Vegetation in all layers 
H_canopy_prof(:,tpq) = H_sun_prof(:,tpq) + H_shade_prof(:,tpq);               %[W/m^2/ground] Total Sensible Heat from Vegetation in all layers 
Rnrad_canopy_prof(:,tpq) = Rnrad_sun_prof(:,tpq) + Rnrad_shade_prof(:,tpq);   %[W/m^2/ground] Net Radiation in the vegetation for all layers  
TR_canopy_prof(:,tpq) = TR_sun_prof(:,tpq) + TR_shade_prof(:,tpq);            %[mm/s] Transpiration for all layers
        
    % Normalized Flux Profiles (ie. per unit LAI)    
    %   Sunlit
An_sun_norm_prof(:,tpq,:) = An_sun;                                         %[umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration sunlit fraction for all layers
LE_sun_norm_prof(:,tpq,:) = LE_sun;                                         %[W/ m^2 leaf / s] Latent Heat sunlit fraction for all layers 
H_sun_norm_prof(:,tpq,:) = H_sun;                                           %[W/ m^2 leaf / s] Sensible Heat sunlit fraction for all layers  
Rnrad_sun_norm_prof(:,tpq) = Rnrad_sun./LAI_sun;                            %[W/ m^2 leaf / s] Net Radiation sunlit fraction for all layers
        
    %   Shaded
An_shade_norm_prof(:,tpq,:) = An_shade;                                     %[umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration shade fraction for all layers
LE_shade_norm_prof(:,tpq,:) = LE_shade;                                     %[W/ m^2 leaf / s] Latent Heat shade fraction for all layers  
H_shade_norm_prof(:,tpq,:) = H_shade;                                       %[W/ m^2 leaf / s] Sensible Heat shade fraction for all layers 
Rnrad_shade_norm_prof(:,tpq) = Rnrad_shade./LAI_shade;                      %[W/ m^2 leaf / s] Net Radiation shade fraction for all layers 
        
    %   Canopy
An_canopy_norm_prof(:,tpq) = (sum(An_sun.*fLAIz,2)).*fsun + (sum(An_shade.*fLAIz,2)).*fshade; % [umol CO2/ m^2 leaf / s] Photosynthesis minus leaf respiration  for all layers
LE_canopy_norm_prof(:,tpq) = (sum(LE_sun.*fLAIz,2)).*fsun + (sum(LE_shade.*fLAIz,2)).*fshade; % [W/ m^2 leaf / s] Latent Heat for all layers
H_canopy_norm_prof(:,tpq) = (sum(H_sun.*fLAIz,2)).*fsun + (sum(H_shade.*fLAIz,2)).*fshade;    % [W/ m^2 leaf / s] Sensible Heat for all layers   
Rnrad_canopy_norm_prof(:,tpq) = (Rnrad_sun.*fsun) + (Rnrad_shade.*fshade);                    % [W/ m^2 leaf / s] Net Radiation for all layers   
            
% CANOPY STATES:        
    % Leaf States
Tl_sun_prof(:,tpq,:) = Tl_sun;                                              %[C] Leaf Temperature Sunlit Fraction all layers and species
Tl_shade_prof(:,tpq,:) = Tl_shade;                                          %[C] Leaf Temperature Shade Fraction all Layers and species 
Tl_sun_Ta_Diff(:,tpq,:) = Tl_sun - repmat(TAz,1,nspecies);                  %[C] Temperature Difference Atmosphere and Leaf Sunlint for all layers and species
Tl_shade_Ta_Diff(:,tpq,:) = Tl_shade - repmat(TAz,1,nspecies);              %[C] Temperature Difference Atmosphere and Leaf Shade for all layers and species
        
psil_sun_prof(:,tpq,:) = psil_sun;                                          %[MPa] Leaf Water Potential sunlit fraction, for all layers and species
psil_shade_prof(:,tpq,:) = psil_shade;                                      %[MPa] Leaf Water Potential Shade fraction, for all layers and species
        
fsvg_sun_prof(:,tpq,:) = fsvg_sun;                                          %[] g Tuzet factor to Ball Berry Model for sunlit fraction, for all layers and species
fsvm_sun_prof(:,tpq,:) = fsvm_sun;                                          %[] m Tuzet factor to Ball Berry Model for sunlit fraction, for all layers and species 
fsvg_shade_prof(:,tpq,:) = fsvg_shade;                                      %[] g Tuzet factor to Ball Berry Model for shade fraction, for all layers and species
fsvm_shade_prof(:,tpq,:) = fsvm_shade;                                      %[] m Tuzet factor to Ball Berry Model for shade fraction, for all layers and species
        
gsv_sun_prof(:,tpq,:) = gsv_sun;                                            %[mol/m^2/s] Stomatal Conductance for sunlit fraction for all layers and species
gsv_shade_prof(:,tpq,:) = gsv_shade;                                        %[mol/m^2/s] Stomatal Conductance for shade fraction for all layers and species
        
Ci_sun_prof(:,tpq,:) = Ci_sun;                                              %[umol/mol] Internal leaf concentration of CO2 for sunlit fraction, for all layers and species
Ci_shade_prof(:,tpq,:) = Ci_shade;                                          %[umol/mol] Internal leaf concentration of CO2 for shade fraction, for all layers and species
        
gbv_sun_prof(:,tpq,:) = gbv_sun;                                            %[mol/m^2/s] Boundary Layer conductance to vapor in sunlit, for all layers and species
gbh_sun_prof(:,tpq,:) = gbh_sun;                                            %[mol/m^2/s] Boundary Layer conductance to heat in sunlit, for all layers and species
gbv_shade_prof(:,tpq,:) = gbv_shade;                                        %[mol/m^2/s] Boundary Layer conductance to vapor in shade, for all layers and species
gbh_shade_prof(:,tpq,:) = gbh_shade;                                        %[mol/m^2/s] Boundary Layer conductance to heat in shade, for all layers and species

LAI_sun_prof(:,tpq) = LAI_sun;                                              %[m^2 leaf / m^2 ground] Leaf Area Index sunlit, for all layers 
LAI_shade_prof(:,tpq) = LAI_shade;                                          %[m^2 leaf / m^2 ground] Leaf Area Index shade, for all layers
        
fsun_prof(:,tpq) = fsun;                                                    %[] Fraction of LAI sunlit
fshade_prof(:,tpq) = fshade;                                                %[] Fraction of LAI shade 
        
LSsunCON_store(1,tpq,:)=LSsunCON;                                           % Convergence of leaf solution, for sunlit, for all species             
LSshaCON_store(1,tpq,:)=LSshaCON;                                           % Convergence of leaf solution, for shade, for all specie 
            
% leaf convergence
cntleafsolution_store(1,tpq,:) = VARIABLES.CANOPY.cntleafsolution;          % Leaf convergence 

% Photosynthetic Biochemistry    
Ph_limit_sun_store(:,tpq,:) = Ph_limit_sun;                                 %[] Type of Photosynthesis 1. Rubisco-limited, 2light-limited, 3. Sucrose Limited
Jc_C3_sun_store(:,tpq,:) = Jc_C3_sun;                                       %[umol/m^2 leaf area/s] Rubisco-Limited Rate for sunlit, for all layers, and species
Jj_C3_sun_store(:,tpq,:) = Jj_C3_sun;                                       %[umol/m^2 leaf area/s] Light-Limited Rate for sunlit, for all layers, and species
Js_C3_sun_store(:,tpq,:) = Js_C3_sun;                                       %[umol/m^2 leaf area/s] Sucrose-Limited Rate for sunlit, for all layers, and species
Jc_C4_sun_store(:,tpq,:) = Jc_C4_sun;                                       %[umol/m^2 leaf area/s] Rubisco-Limited Rate for sunlit, for all layers, and species
Jj_C4_sun_store(:,tpq,:) = Jj_C4_sun;                                       %[umol/m^2 leaf area/s] Ligth-Limited Rate for sunlit, for all layers, and species   
Js_C4_sun_store(:,tpq,:) = Js_C4_sun;                                       %[umol/m^2 leaf area/s] Sucrose-Limited Rate for sunlit, for all layers, and species
        
Ph_limit_shade_store(:,tpq,:) = Ph_limit_shade;                             %[] Type of Photosynthesis 1. Rubisco-limited, 2light-limited, 3. Sucrose Limited 
Jc_C3_shade_store(:,tpq,:) = Jc_C3_shade;                                   %[umol/m^2 leaf area/s] Rubisco-Limited Rate for shade, for all layers, and species                                 
Jj_C3_shade_store(:,tpq,:) = Jj_C3_shade;                                   %[umol/m^2 leaf area/s] Ligth-Limited Rate for shade, for all layers, and species 
Js_C3_shade_store(:,tpq,:) = Js_C3_shade;                                   %[umol/m^2 leaf area/s] Sucrose-Limited Rate for shade, for all layers, and species 
Jc_C4_shade_store(:,tpq,:) = Jc_C4_shade;                                   %[umol/m^2 leaf area/s] Rubisco-Limited Rate for shade, for all layers, and species 
Jj_C4_shade_store(:,tpq,:) = Jj_C4_shade;                                   %[umol/m^2 leaf area/s] Ligth-Limited Rate for shade, for all layers, and species 
Js_C4_shade_store(:,tpq,:) = Js_C4_shade;                                   %[umol/m^2 leaf area/s] Sucrose-Limited Rate for shade, for all layers, and species        
                    
dryfrac_prof(:,tpq) = dryfrac;                                              %[] Fraction of LAI that is dry
wetfrac_prof(:,tpq) = wetfrac;                                              %[] Fraction of LAI that is wet
                
ppt_ground_store(:,tpq,:) = ppt_ground;                                       %[mm] Precipitation Reaching the ground
    
    % Mean Canopy Profiles    
Ci_canopy_prof(:,tpq) = (nansum(Ci_sun.*fLAIz,2)).*fsun + (nansum(Ci_shade.*fLAIz,2)).*fshade;       %[umol/mol] Internal leaf concentration of CO2, for all layers
Tl_canopy_prof(:,tpq) = (nansum(Tl_sun.*fLAIz,2)).*fsun + (nansum(Tl_shade.*fLAIz,2)).*fshade;       %[C] Leaf Temperature, for all layers 
gsv_canopy_prof(:,tpq) = (nansum(gsv_sun.*fLAIz,2)).*fsun + (nansum(gsv_shade.*fLAIz,2)).*fshade;    %[mol/m^2/s] Stomatal Conductance for all layers  
psil_canopy_prof(:,tpq) = (nansum(psil_sun.*fLAIz,2)).*fsun + (nansum(psil_shade.*fLAIz,2)).*fshade; %[MPa] Leaf Water Potential for all layers 
fsvg_canopy_prof(:,tpq) = (nansum(fsvg_sun.*fLAIz,2)).*fsun + (nansum(fsvg_shade.*fLAIz,2)).*fshade; %[] g Tuzet factor to Ball Berry Model for all layers and species
fsvm_canopy_prof(:,tpq) = (nansum(fsvm_sun.*fLAIz,2)).*fsun + (nansum(fsvm_shade.*fLAIz,2)).*fshade; %[] m Tuzet factor to Ball Berry Model for all layers and species 
gbv_canopy_prof(:,tpq) = (nansum(gbv_sun.*fLAIz,2)).*fsun + (nansum(gbv_shade.*fLAIz,2)).*fshade;    %[mol/m^2/s] Boundary Layer conductance to vapor for all layers 
gbh_canopy_prof(:,tpq) = (nansum(gbh_sun.*fLAIz,2)).*fsun + (nansum(gbh_shade.*fLAIz,2)).*fshade;    %[mol/m^2/s] Boundary Layer conductance to heat for all layers 
                
    % Mean Canopy States
Ci_mean(1,tpq) = sum(Ci_canopy_prof(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);                         %[umol/mol] Internal leaf concentration of CO2, mean canopy
Tl_mean(1,tpq) = sum(Tl_canopy_prof(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);                         %[C] Leaf Temperature of CO2, mean canopy
gsv_mean(1,tpq) = sum(gsv_canopy_prof(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);                       %[mol/m^2/s] Stomatal Conductance, mean canopy  
psil_mean(1,tpq) = sum(psil_canopy_prof(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);                     %[MPa] Leaf Water Potential for all layers
fsvg_mean(1,tpq) = sum(fsvg_canopy_prof(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);                     %[] g Tuzet factor to Ball Berry Model for all layers and specie
fsvm_mean(1,tpq) = sum(fsvm_canopy_prof(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);                     %[] m Tuzet factor to Ball Berry Model for all layers and specie


    % Canopy Microenvironment
CAz_prof(:,tpq,:) = VARIABLES.CANOPY.CAz;                                                      %[umol/mol] Atmospheric Concentration of CO2, for all layers
TAz_prof(:,tpq,:) = VARIABLES.CANOPY.TAz;                                                      %[C] Air Temperature, for all layers
EAz_prof(:,tpq,:) = VARIABLES.CANOPY.EAz;                                                      %[kPa] Vapor Pressure Air, for all layers  
Uz_prof(:,tpq) = VARIABLES.CANOPY.Uz;                                                        %[m/s] Wind Velocity
                
    % Canopy H2O Storage
Sh2o_canopy_prof(:,tpq) = Sh2o_prof;                                        %[mm] Intercepted water in canopy, for all layers 
Sh2o_canopy(1,tpq) = Sh2o_can;                                              %[mm] Total Intercepted water in canopy.
VARIABLES.CANOPY.Sh2o_can_prev = Sh2o_can;                                 %[mm] Total Intercepted water in canopy in previous time step 
        
    % Condensation
Ch2o_canopy_prof(:,tpq) = Ch2o_prof;                                        %[mm] Condensation water in the canopy, for all layers  
Ch2o_canopy(1,tpq) = Ch2o_can;                                              %[mm] Total Condensation water in the canopy.  
        
    % Evaporation
Evap_canopy_prof(:,tpq) = Evap_prof;                                        %[mm] Evaporation water from canopy 
Evap_canopy(1,tpq) = Evap_can;                                              %[mm] Total Evaporation water from canopy.


% SOIL VARIABLES
volliq_store(:,tpq,:) = VARIABLES.SOIL.volliq;                                               %[] volumetric water content  
krad_store(:,tpq,:) = VARIABLES.ROOT.krad;                                                 %[/s] radial conductivity of the root
kax_store(:,tpq,:) = VARIABLES.ROOT.kax;                                                   %[mm/s] axial conductivity of the root    
hk_store(:,tpq) = hk;                                                       %[mm/s] Soil Hydraulic Conductivity at Node Layers
smp_store(:,tpq,:) = VARIABLES.SOIL.smp;                                                     %[mm of H2O] Soil Matric Potential 
rpp_store(:,tpq,:) = VARIABLES.ROOT.rpp;                                                   %[mm of H2O] Root Water Potential 
mberrormm_store(1,tpq)=mberrormm;                                           %[mm/dtime] Mass balance error.

    
smpMPA_store(:,tpq,:) = VARIABLES.SOIL.smp*CONSTANTS.mmH2OtoMPa;                             %[MPa] Soil Matric Potential for all layers
rppMPA_store(:,tpq,:) = VARIABLES.ROOT.rpp*CONSTANTS.mmH2OtoMPa;                           %[MPa] Root Water Potential for all layers and species
      
   % smp_weight_store(tpq) = VARIABLES.SOIL.smp_wgt*CONSTANTS.mmH2OtoMP                     %[MPA] Weighted mean smp over root uptake profile   
    rpp_weight_store(1,tpq,:) = VARIABLES.ROOT.rpp_wgt*CONSTANTS.mmH2OtoMPa;                %[Mpa] Weighted mean Root weight to use in Canopy Model      
%    volliq_weight_store(tpq) = sum(VARIABLES.ROOT.volliq .* VERTSTRUC.rootfr);
    
qlayer_store(:,tpq) = VARIABLES.SOIL.qlayer;                                %[mm/s] Fluxes of water between layers in the soil
bot_drainage_store(1,tpq) = VARIABLES.SOIL.qlayer(end);                                    %[mm/s] Flux of water at bottom of soil column
hor_drainage_store(1,tpq) = VARIABLES.SOIL.hor_drainage;                    %[mm/s] Horizontal Drainage. Only happens when superaturation occurs.   
hor_drainage_lay_store(:,tpq) = VARIABLES.SOIL.hor_drainage_lay;            %[mm/s] Horizontal Drainage for all layers. Only happens when superaturation occurs.   
    
    
if (SWITCHES.soilheat_on)
        Ts_store(:,tpq,:) = VARIABLES.SOIL.Ts;                                %[C] Temperature in the soil. For all layers  
        Gnew_store(:,tpq,:) = VARIABLES.SOIL.Gnew;                            %[C] Soil Heat Flux. For all layers in the soil. 
        cpv_store(:,tpq,:) = VARIABLES.SOIL.cpv;                              %[J / m^3 / K] Soil Volumetric Heat Capacity  
        TKsoil_store(:,tpq,:) = VARIABLES.SOIL.TKsoil;                        %[W / m / K] Thermal Conductivity Soil 
        TKsol_store(:,tpq,:) = VARIABLES.SOIL.TKsol;                          %[W / m / K] Thermal Conductivity of Mineral Solids   
end
    
if (SWITCHES.ns)
        wuptake_store(:,tpq) = layeruptake;                                 %[mm/s] Fluxes of plant water uptake at each layer. Shows the total from all the species 
        wuptake_all_store(:,tpq,:) =layeruptake_all;                        %[mm/s] Fluxes of plant water uptake at each layer, and for all species. 
end


runoff_store(1,tpq,:) = VARIABLES.SOIL.runoff;                                %[mm] Flux of water lost as runoff 
qback_store(1,tpq,:) = VARIABLES.SOIL.qback;                                  %[mm/s] Flux of water returned to litter or surface by sypersaturation of first layer
qadflux_store(1,tpq,:) = VARIABLES.SOIL.qadflux;                              %[mm/s]  Final Flux that Infiltrates in the soil
qback_lit_store(1,tpq,:) = VARIABLES.SOIL.qback_lit;                          %[mm/s] Real flux that is going back to snow-litter pack  
    
Gs_store(1,tpq,:) = VARIABLES.SOIL.Gs;                                        %Ground Heat Flux into the Soil from the Soil- Snow-Litter Pack  Boundary [W/m2]    
Gsl_store(1,tpq,:) = VARIABLES.SOIL.Gsl;                                      %Ground Heat Flux into the Soil- Snow Litter Pack Boundary from The Snow-Litter Pack [W/m2] 
LE_soli_store(1,tpq,:) = VARIABLES.SOIL.LEs;                                  %[W/m2] Latent Heat from the soil only  
LE_liat_store(1,tpq,:) = VARIABLES.SOIL.LEl;                                  %[W/m2] Latent Heat from snow litter pack  
psili_MPa_store(1,tpq,:) = VARIABLES.SOIL.psil_MPa;                           %[Mpa] Soil Water Potential in the snow-Litter Pack [converted from m to mm]
psili_store(1,tpq,:) = VARIABLES.SOIL.psili;                                  %[mm] Soil Water Potential in the snow-Litter Pack 
Tli_store(1,tpq,:) = VARIABLES.SOIL.Tli;                                      %[C] Temperature Snow-Litter pack  
Tsl_store(1,tpq,:) = VARIABLES.SOIL.Tsl;                                      %[C] Tempeature in the Soil-Litter Interface 
    % rhosn
rhosn_store(1,tpq,:) = VARIABLES.SOIL.rhosn;                                  %[kg/m3] Snow Density
    % w
wliqsl_store(1,tpq,:) = VARIABLES.SOIL.wliqsl;                                %[kg / m^2] Liquid Water Density (per unit area) 
wicesl_store(1,tpq,:) = VARIABLES.SOIL.wicesl;                                %[kg / m^2] Ice Water Density (per unit area) 
deltawice_store(1,tpq,:) = VARIABLES.SOIL.deltawice;                          %[] Delta of ice needed for change of phase in the snow-litter pack
dH_store(1,tpq,:) = VARIABLES.SOIL.dH;                                        %[W/m2] Delta of Energy Melting/Fusion
dS_store(1,tpq,:) = VARIABLES.SOIL.dS;                                        %[W/m2] Change of Internal Energy in Snow-Litter Pack in Time Step
    
wsn_store(1,tpq,:) = VARIABLES.SOIL.wsn;                                      %[kg / m^2] Snow Density (per unit area) 
    % z
zicesl_store(1,tpq,:) = VARIABLES.SOIL.zicesl;                                %[mm] Ice water depth in Snow-litter pack 
zliqsl_store(1,tpq,:) = VARIABLES.SOIL.zliqsl;                                %[mm] Liquid water depth  in Snow-litter pack
%VARIABLES.SOIL.zliqsl_prev = zliqsl_store(max(tpq,1));                     %[mm] Ice water depth in Snow-litter pack in previous time step     
%VARIABLES.SOIL.zicesl_prev = zicesl_store(max(tpq,1));                     %[mm] Liquid water depth in Snow-litter pack in previous time step
    
zsn_store(1,tpq,:) = VARIABLES.SOIL.zsn;                                      %[mm] Depth of snow 
    % vols
volliqli_store(1,tpq,:) = VARIABLES.SOIL.volliqli ;                           %[] Volumetric Water Content in Litter Pack   
voliceli_store(1,tpq,:) = VARIABLES.SOIL.voliceli;                            %[] Volumetric Ice Content in Litter Pack     
volliqsn_store(1,tpq,:) = VARIABLES.SOIL.volliqsn;                            %[] Volumetric Liquid Water Content in Snow    
volicesn_store(1,tpq,:) = VARIABLES.SOIL.volicesn;                            %[] Volumetric Ice Content in Snow        
voltotsn_store(1,tpq,:) = VARIABLES.SOIL.voltotsn;                            %[] Total Volumetric Ice and Liquid Water Content in Snow      
voltotli_store(1,tpq,:) = VARIABLES.SOIL.voltotli;                            %[] Total Volumetric Ice and Liquid Water Content in Litter         
    %Other variables 
CRm_store(1,tpq,:) = VARIABLES.SOIL.CRm;                                      %[] Snow Compaction Factor 
CRo_store(1,tpq,:) = VARIABLES.SOIL.CRo;                                      %[] Snow Compaction Factor 
drhosn_store(1,tpq,:) = VARIABLES.SOIL.drhosn;                                %[kg / m^3] Change in snow density        
qinflL_store(1,tpq,:) = VARIABLES.SOIL.qinflL;                                %[mm/s] Water that infiltrates into snow 
net_qinflL_store(1,tpq,:) = VARIABLES.SOIL.net_qinflL;                        %[mm/s] Net water that infiltrates into snow excluding drainage in soil
qinfl_store(1,tpq,:) = VARIABLES.SOIL.qinfl;                                  %[mm/s] Infiltration into soil 
drainlitter_store(1,tpq,:) = VARIABLES.SOIL.drainlitter;                      %[mm] Water that drains into soil from litter 
pptrate_ground_store(1,tpq,:) = VARIABLES.SOIL.pptrate_ground;                %[mm/s] Precipitation reaching the ground 
Esoil_store(1,tpq,:) = VARIABLES.SOIL.Esoil;                                  %[mm/s] Evaporation from soil  
Esl_store(1,tpq,:) = VARIABLES.SOIL.Esl;                                      %[mm/s]Evaporation or Sublimation from Snow-Litter pack  
snow_tcount_store(1,tpq,:) = VARIABLES.SOIL.snow_tcount;                      %[] Number of time steps with snow to account for snow age 


%% Rohit start
% CROP GROWTH

LeafCMass_can_store(1,tpq,:) = CROP_GROWTH.LeafMass_can;
StemCMass_can_store(1,tpq,:) = CROP_GROWTH.StemMass_can;
RootCMass_can_store(1,tpq,:) = CROP_GROWTH.RootMass_can;
GrainCMass_can_store(1,tpq,:) = CROP_GROWTH.GrainMass_can;

LAIsim_can_store(1,tpq,:) = CROP_GROWTH.LAIsim_can;

TotalCMass_can_store(1,tpq,:) = CROP_GROWTH.LeafMass_can + CROP_GROWTH.StemMass_can + CROP_GROWTH.RootMass_can + CROP_GROWTH.GrainMass_can;
AboveCMass_can_store(1,tpq,:) = CROP_GROWTH.LeafMass_can + CROP_GROWTH.StemMass_can + CROP_GROWTH.GrainMass_can;


LeafBiomass_can_store(1,tpq,:) = CROP_GROWTH.LeafBiomass;
StemBiomass_can_store(1,tpq,:) = CROP_GROWTH.StemBiomass;
RootBiomass_can_store(1,tpq,:) = CROP_GROWTH.RootBiomass;
GrainBiomass_can_store(1,tpq,:) = CROP_GROWTH.GrainBiomass;

TotalBioMass_can_store(1,tpq,:) = CROP_GROWTH.TotalBiomass_can;
AboveBioMass_can_store(1,tpq,:) = CROP_GROWTH.AboveBiomass_can;

GrainYield_store(1,tpq,:) = CROP_GROWTH.GrainYield;

DeathRate_store(1,tpq,:) = VARIABLES.CANOPY.DeathRate;
NDemand_store(1,tpq,:) = VARIABLES.CANOPY.N_Demand;

Nreq_store(1,tpq,:) = VARIABLES.CANOPY.Net_Nreq;
Fert_req_store(1,tpq,:) = VARIABLES.CANOPY.Fert_req;
Fert_doy_store(1,tpq,:) = VARIABLES.CANOPY.N_Fert_DOY;

K_ec_store(1,tpq,:)= VARIABLES.NCan.K_ec;
K_rub_leaf_store(1,tpq,:)= VARIABLES.NCan.K_rub_leaf;
K_lhc_leaf_store(1,tpq,:)= VARIABLES.NCan.K_lhc_leaf;
fn_store(1,tpq,:)=VARIABLES.NCan.fn;

% Irrigation

IW_store(1,tpq,:) = VARIABLES.SOIL.IW;

VARIABLES.CANOPY.psil_mean = sum(psil_canopy_prof(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);
VARIABLES.CANOPY.gsv_mean = sum(gsv_canopy_prof(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);
%
% Tl_canopy_prof_up(:,tpq) = (nansum(Tl_sun_up.*fLAIz,2)).*fsun + (nansum(Tl_shade_up.*fLAIz,2)).*fshade;
%
% Tl_canopy_prof_lo(:,tpq) = (nansum(Tl_sun_lo.*fLAIz,2)).*fsun + (nansum(Tl_shade_lo.*fLAIz,2)).*fshade;
%
%
Tl_mean_up_store(1,tpq,:) = VARIABLES.CANOPY.Tl_mean_up;
Tl_mean_lo_store(1,tpq,:) = VARIABLES.CANOPY.Tl_mean_lo;

% VARIABLES.CANOPY.Tl_mean_up = sum(Tl_canopy_prof_up(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);
%
% VARIABLES.CANOPY.Tl_mean_lo = sum(Tl_canopy_prof_lo(vinds,tpq).*LADnorm(vinds))/sum(LADnorm);

CWSI_store(1,tpq,:) = CROP_GROWTH.CWSI;

% Nitrogen Canopy

% Vcmax_store(1,tpq) = VARIABLES.CANOPY.Vcmax;
Vcmax_vert_store(:,tpq, :) = VARIABLES.NCan.Vcmax_vert_x;
Chl_la_store(1,tpq,:) = VARIABLES.NCan.Chl;
N_can_store(1,tpq,:) = VARIABLES.NCan.N_can;
N_leaf_la_store(:,tpq,:) = VARIABLES.NCan.Nleaf;
N_rub_la_store(:,tpq,:) = VARIABLES.NCan.Nrub;

N_uptake_can_store(1,tpq,:) = VARIABLES.NCan.N_uptake_can;

Nl_ncan_ga_store(1,tpq,:) = VARIABLES.NCan.Nl;
Nr_ncan_ga_store(1,tpq,:) = VARIABLES.NCan.Nr;
Ns_ncan_ga_store(1,tpq,:) = VARIABLES.NCan.Ns;
Ng_ncan_ga_store(1,tpq,:) = VARIABLES.NCan.Ng;

%% Rohit end

    
% MASS BALANCE 

if tpq>1    
    %  Mass Balance Error 
    MBerror_soil(1,tpq,:) = VARIABLES.MB.mbsoil/86400;                        % [mm/s]
    MBerror_littersoil(1,tpq,:) = VARIABLES.MB.mblittersoil/86400;            % [mm/s]
    MBerror_mbcan(1,tpq,:) = VARIABLES.MB.mbcan/86400;                        % [mm/s]
    MBerror_mbcanlittersoil(1,tpq,:) = VARIABLES.MB.mbcanlittersoil/86400;    % [mm/s]        
%     %  Save Mass Balance Structure
%     MB_store.PPTf(1,tpq) = VARIABLES.MB.PPTf;                               %[mm/s] Total Rainfall  
%     MB_store.PPTgr(1,tpq) = VARIABLES.MB.PPTgr;                             %[mm/s] Total Rainfall Reaching Ground  
%     MB_store.INFsoil(1,tpq) = VARIABLES.MB.INFsoil;                         %[mm/s] Infiltration into Soil, Drainage From Litter 
%     MB_store.EVsoil(1,tpq) = VARIABLES.MB.EVsoil;                           %[mm/s] Evaporation from soil 
%     MB_store.EVcanf(1,tpq) = VARIABLES.MB.EVcanf;                           %[mm/s] Evaporation from canopy 
%     MB_store.Bdrainage(1,tpq) = VARIABLES.MB.Bdrainage;                     %[mm/s] Bottom Drainage from soil 
%     MB_store.TR(1,tpq,:) = VARIABLES.MB.TR;                                 %[mm/s] Transpiration all species 
%     MB_store.runoff(1,tpq) = VARIABLES.MB.runoff;                           %[mm/s] Transpiration all species 
%     MB_store.COcanf(1,tpq) = VARIABLES.MB.COcanf;                           %[mm/s] Total canopy condensation 
% 
%     MB_store.dws_f(1,tpq) = VARIABLES.MB.dws_f;                             %[mm/s] Rate of change in soil water content 
%     MB_store.dwc_f(1,tpq) = VARIABLES.MB.dwc_f;                             %[mm/s] Rate of change in canopy water content 
%     MB_store.dwl_f(1,tpq) = VARIABLES.MB.dwl_f;                             %[mm/s] Rate of change in litter water content 
% 
%     MB_store.mbsoil(1,tpq) = VARIABLES.MB.mbsoil;                           %[mm/day] Mass balance error soil 
%     MB_store.mbcan(1,tpq) = VARIABLES.MB.mbcan;                             %[mm/day] Mass balance error canopy     
%     MB_store.mblittersoil(1,tpq) = VARIABLES.MB.mblittersoil;               %[mm/day] Mass balance error litter-soil 
%     MB_store.mbcanlittersoil(1,tpq) = VARIABLES.MB.mbcanlittersoil;         %[mm/day] Mass balance error litter-soil-canopy
end



if (SWITCHES.entropy_on)
        % CANOPY
        % in
  SScandir_in_all_lay_store(:,tpq,:) = SSresults.Total.SScandir_in_all_lay; %[W/m^2/K] Entropy in incoming direct SW in canopy for all layers and species
  SScandir_in_all_store(1,tpq,:) = SSresults.Total.SScandir_in_all;         %[W/m^2/K] Entropy in incoming direct SW in canopy for all layers  
  SScandir_in_tot_store(1,tpq) = SSresults.Total.SScandir_in_tot;           %[W/m^2/K] Entropy in incoming direct SW in canopy   

  SScandif_in_all_lay_store(:,tpq,:) = SSresults.Total.SScandif_in_all_lay; %[W/m^2/K] Entropy in incoming diffuse SW in canopy for all layers and species
  SScandif_in_all_store(1,tpq,:) = SSresults.Total.SScandif_in_all;         %[W/m^2/K] Entropy in incoming diffuse SW in canopy for all layers 
  SScandif_in_tot_store(1,tpq) = SSresults.Total.SScandif_in_tot;           %[W/m^2/K] Entropy in incoming diffuse SW in canopy  

  SScanLW_in_all_lay_store(:,tpq,:) = SSresults.Total.SScanLW_in_all_lay;   %[W/m^2/K] Entropy in incoming LW in canopy for all layers and species
  SScanLW_in_all_store(1,tpq,:) = SSresults.Total.SScanLW_in_all;           %[W/m^2/K] Entropy in incoming LW in canopy for all layers 
  SScanLW_in_tot_store(1,tpq) = SSresults.Total.SScanLW_in_tot;             %[W/m^2/K] Entropy in incoming LW in canopy 

  SScan_in_all_lay_store(:,tpq,:) = SSresults.Total.SScan_in_all_lay;       %[W/m^2/K] Entropy in total energy fluxes leaving canopy for all layers and species
  SScan_in_all_store(1,tpq,:) = SSresults.Total.SScan_in_all;               %[W/m^2/K] Entropy in total energy fluxes leaving canopy for all layers  
  SScan_in_store(1,tpq) = SSresults.Total.SScan_in;                         %[W/m^2/K] Entropy in total energy fluxes leaving canopy   

        % out
  SScandir_out_all_lay_store(:,tpq,:) = SSresults.Total.SScandir_out_all_lay;%[W/m^2/K] Entropy in outgoing direct SW in canopy for all layers and species
  SScandir_out_all_store(1,tpq,:) = SSresults.Total.SScandir_out_all;       %[W/m^2/K] Entropy in outgoing direct SW in canopy for all layers  
  SScandir_out_tot_store(1,tpq) = SSresults.Total.SScandir_out_tot;         %[W/m^2/K] Entropy in outgoing direct SW in canopy  

  SScandif_out_all_lay_store(:,tpq,:) = SSresults.Total.SScandif_out_all_lay;%[W/m^2/K] Entropy in outgoing diffuse SW in canopy for all layers and species
  SScandif_out_all_store(1,tpq,:) = SSresults.Total.SScandif_out_all;       %[W/m^2/K] Entropy in outgoing diffuse SW in canopy for all layers  
  SScandif_out_tot_store(1,tpq) = SSresults.Total.SScandif_out_tot ;        %[W/m^2/K] Entropy in outgoing diffuse SW in canopy  

  SScanLW_out_all_lay_store(:,tpq,:) = SSresults.Total.SScanLW_out_all_lay; %[W/m^2/K] Entropy in outgoing diffuse LW in canopy for all layers and species
  SScanLW_out_all_store(1,tpq,:) = SSresults.Total.SScanLW_out_all;         %[W/m^2/K] Entropy in outgoing diffuse LW in canopy for all layers  
  SScanLW_out_tot_store(1,tpq) = SSresults.Total.SScanLW_out_tot;           %[W/m^2/K] Entropy in outgoing diffuse LW in canopy  

  SScanLE_out_all_lay_store(:,tpq,:) = SSresults.Total.SScanLE_out_all_lay; %[W/m^2/K] Entropy associated with LE in canopy for all layers and species
  SScanLE_out_all_store(1,tpq,:) = SSresults.Total.SScanLE_out_all;         %[W/m^2/K] Entropy associated with LE in canopy for all layers 
  SScanLE_out_tot_store(1,tpq) = SSresults.Total.SScanLE_out_tot;           %[W/m^2/K] Entropy associated with LE in canopy  

  SScanH_out_all_lay_store(:,tpq,:) = SSresults.Total.SScanH_out_all_lay;   %[W/m^2/K] Entropy associated with H in canopy for all layers and species 
  SScanH_out_all_store(1,tpq,:) = SSresults.Total.SScanH_out_all;           %[W/m^2/K] Entropy associated with H in canopy for all layers  
  SScanH_out_tot_store(1,tpq) = SSresults.Total.SScanH_out_tot;             %[W/m^2/K] Entropy associated with H in canopy 

        % total out
  SScan_out_all_lay_store(:,tpq,:) = SSresults.Total.SScan_out_all_lay;     %[W/m^2/K] Entropy in total energy fluxes leaving canopy for all layers and species  
  SScan_out_all_store(1,tpq,:) = SSresults.Total.SScan_out_all;             %[W/m^2/K] Entropy in total energy fluxes leaving canopy for all layers  
  SScan_out_store(1,tpq) = SSresults.Total.SScan_out;                       %[W/m^2/K] Entropy in total energy fluxes leaving canopy  

        % TOTAL CANOPY
  SScan_all_lay_store(:,tpq,:) = SSresults.Total.SScan_all_lay;             %[W/m^2/K] Net Entropy Budget in Canopy for all layers and species  
  SScan_all_store(1,tpq,:) = SSresults.Total.SScan_all;                     %[W/m^2/K] Net Entropy Budget in Canopy for all layers  
  SScan_tot_store(1,tpq) = SSresults.Total.SScan_tot;                       %[W/m^2/K] Net Entropy Budget in Canopy   
        
        % PHOTOSYNTHESIS
  Epho_dif_store(1,tpq) = SSresults.Epho_dif;                               %[W/m^2] Energy associated with diffuse SW used in photosynthesis  
  Epho_dir_store(1,tpq) = SSresults.Epho_dir;                               %[W/m^2] Energy associated with direct SW used in photosynthesis 
  Epho_tot_store(1,tpq) = SSresults.Epho_tot;                               %[W/m^2] Entropy associated with SW used in photosynthesis 
  SSphodir_in_store(1,tpq) = SSresults.SSphodir_in;                         %[W/m^2/K] Entropy associated with direct SW used in photosynthesis                            
  SSphodif_in_store(1,tpq) = SSresults.SSphodif_in;                         %[W/m^2/K] Entropy associated with diffuse SW used in photosynthesis 
                
        % SOIL
        % in
  SSsoildir_in_tot_store(1,tpq) = SSresults.Total.SSsoildir_in_tot;         %[W/m^2/K] Entropy associated with direct SW reaching soil 
  SSsoildif_in_tot_store(1,tpq) = SSresults.Total.SSsoildif_in_tot;         %[W/m^2/K] Entropy associated with diffuse SW reaching soil 
  SSsoilLW_in_tot_store(1,tpq) = SSresults.Total.SSsoilLW_in_tot;           %[W/m^2/K] Entropy associated with LW reaching soil 
  SSsoil_in_store(1,tpq) = SSresults.Total.SSsoil_in;                       %[W/m^2/K] Entropy in total energy fluxes reaching soil    
        
        % out 
  SSsoildir_out_tot_store(1,tpq) = SSresults.Total.SSsoildir_out_tot;       %[W/m^2/K] Entropy associated with direct SW leaving soil       
  SSsoildif_out_tot_store(1,tpq) = SSresults.Total.SSsoildif_out_tot;       %[W/m^2/K] Entropy associated with diffuse SW leaving soil 
  SSsoilLW_out_tot_store(1,tpq) = SSresults.Total.SSsoilLW_out_tot;         %[W/m^2/K] Entropy associated with LW leaving soil 
  SSsoilLE_out_tot_store(1,tpq) = SSresults.Total.SSsoilLE_out_tot;         %[W/m^2/K] Entropy associated with LE dissipated from soil  
  SSsoilH_out_tot_store(1,tpq) = SSresults.Total.SSsoilH_out_tot;           %[W/m^2/K] Entropy associated with H dissipated from soil     
  SSsoil_out_store(1,tpq) = SSresults.Total.SSsoil_out;                     %[W/m^2/K] Entropy in total energy fluxes leaving soil 
        
        % dH CANOPY
  SSdHcan_store(1,tpq) = SSresults.SSdHcan_tot;                             %[W/m^2/K] Entropy associated with heat stored in canopy  
                               
        % G
  SSsoilG_store(1,tpq) = SSresults.SSsoilG;                                 %[W/m^2/K] Entropy associated with ground heat flux
        
        % dH dS SOIL
  SSsoildH_store(1,tpq) = SSresults.SSsoildH;                               %[W/m^2/K] Entropy associated with dH in snow (heat for melting/freezing)  
  SSsoildS_store(1,tpq) = SSresults.SSsoildS;                               %[W/m^2/K] Entropy associated with dS in snow (heat stored in snow-litter pack) 
        
        
        % TOTAL SOIL
  SSsoil_tot_store(1,tpq) = SSresults.Total.SSsoil_tot;                     %[W/m^2/K] Net Entropy Budget in Soils  

        % TOTAL ECOSYSTEM
  SSeco_tot_store(1,tpq) = SSresults.Total.SSeco_tot;                       %[W/m^2/K] Net Entropy Budget in Ecosystem  
        
        % TOTALS BY TYPE OF ENERGY
  % SW direct
  SSeco_totSWdir_store(1,tpq) = (SSresults.Total.SScandir_in_tot - SSresults.Total.SScandir_out_tot) + (SSresults.Total.SSsoildir_in_tot - SSresults.Total.SSsoildir_out_tot);
  % SW diffuse
  SSeco_totSWdif_store(1,tpq) = (SSresults.Total.SScandif_in_tot - SSresults.Total.SScandif_out_tot) + (SSresults.Total.SSsoildif_in_tot - SSresults.Total.SSsoildif_out_tot);
  % LW   
  SSeco_totLW_store(1,tpq) = (SSresults.Total.SScanLW_in_tot - SSresults.Total.SScanLW_out_tot) + (SSresults.Total.SSsoilLW_in_tot - SSresults.Total.SSsoilLW_out_tot);
  % H  
  SSeco_totH_store(1,tpq) = -(SSresults.Total.SScanH_out_tot + SSresults.Total.SSsoilH_out_tot);
  % LE  
  SSeco_totLE_store(1,tpq) = -(SSresults.Total.SScanLE_out_tot + SSresults.Total.SSsoilLE_out_tot);
  % TOT (xxx to test: should be equal to previous total xxx)  
  SSeco_tot_store_test(1,tpq) = SSeco_totSWdir_store(1,tpq) + SSeco_totSWdif_store(1,tpq) + SSeco_totLW_store(1,tpq) + SSeco_totH_store(1,tpq) + SSeco_totLE_store(1,tpq);
        
        
        % NET ENERGY
  SWdir_in_store(1,tpq) = SSresults.SWdir_in;                               % [W/m^2] Net (All Ecosystem) Energy incoming SW direct  
  SWdir_out_store(1,tpq) = SSresults.SWdir_out;                             % [W/m^2] Net (All Ecosystem) Energy outgoing SW direct 
  SWdif_in_store(1,tpq) = SSresults.SWdif_in;                               % [W/m^2] Net (All Ecosystem) Energy incoming SW diffuse 
  SWdif_out_store(1,tpq) = SSresults.SWdif_out;                             % [W/m^2] Net (All Ecosystem) Energy outgoing SW diffuse 
  LWin_net_store(1,tpq) = SSresults.LWin_net;                               % [W/m^2] Net (All Ecosystem) Energy incoming LW  
  LWout_net_store(1,tpq) = SSresults.LWout_net;                             % [W/m^2] Net (All Ecosystem) Energy outgoing LW  
  LWemi_net_store(1,tpq) = SSresults.LWemi_net;                             % [W/m^2] Net (All Ecosystem) Energy emitted LW  
  LE_net_store(1,tpq) = SSresults.LE_net;                                   % [W/m^2] Net (All Ecosystem) Energy dissipated LE  
  H_net_store(1,tpq) = SSresults.H_net;                                     % [W/m^2] Net (All Ecosystem) Energy dissipated LW  
  G_net_store(1,tpq) = SSresults.G_net;                                     % [W/m^2] Net (All Ecosystem) Energy ground heat flux 
  dHcan_net_store(1,tpq) = SSresults.dHcan_net;                             % [W/m^2] Net (All Ecosystem) Energy storage in canopy 
                
        % TOTAL ENTROPY USING A NET TEMPERATURE
  SSdir_net_in_store(1,tpq) = SSresults.SSnetdir_in;                        % [W/m^2/K] Net (All Ecosystem) Entropy, due to incoming SW direct 
  SSdir_net_out_store(1,tpq) = SSresults.SSnetdir_out;                      % [W/m^2/K] Net (All Ecosystem) Entropy, due to outgoing SW direct 
  SSdif_net_in_store(1,tpq) = SSresults.SSnetdif_in;                        % [W/m^2/K] Net (All Ecosystem) Entropy, due to incoming SW diffuse 
  SSdif_net_out_store(1,tpq) = SSresults.SSnetdif_out;                      % [W/m^2/K] Net (All Ecosystem) Entropy, due to outgoing SW diffuse 
  SSSW_net_inX_store(1,tpq) = SSresults.SSnetSW_Xin;                        % [] Net (All Ecosystem) X_SW function for computation of entropy in incoming SW radiation 
  SSSW_net_outX_store(1,tpq) = SSresults.SSnetSW_Xout;                      % [] Net (All Ecosystem) X_SW function for computation of entropy in outgoing SW radiation 

  SSLW_net_in_store(1,tpq) = SSresults.SSnetLW_in;                          % [W/m^2/K] Net (All Ecosystem) Entropy, due to incoming LW 
  SSLW_net_inX_store(1,tpq) = SSresults.SSnetLW_Xin;                        % [] Net (All Ecosystem) X_LW function for computation of entropy in incoming LW radiation 
  SSLW_net_out_store(1,tpq) = SSresults.SSnetLW_out;                        % [W/m^2/K] Net (All Ecosystem) Entropy, due to outgoing LW  
  SSLW_net_outX_store(1,tpq) = SSresults.SSnetLW_Xout;                      % [] Net (All Ecosystem) X_LW function for computation of entropy in outgoing LW radiation
        
                
  SSLE_net_store(1,tpq) = SSresults.SSLE_net;                               % [W/m^2/K] Net (All Ecosystem) Entropy, due to LE  
  SSH_net_store(1,tpq) = SSresults.SSH_net;                                 % [W/m^2/K] Net (All Ecosystem) Entropy, due to H  
  SSG_net_store(1,tpq) = SSresults.SSG_net;                                 % [W/m^2/K] Net (All Ecosystem) Entropy, due to G  
  SSdHcan_net_store(1,tpq) = SSresults.SSdHcan_net;                         % [W/m^2/K] Net (All Ecosystem) Entropy, due to storage of heat in canopy  
                
  SSeco_net_in_store(1,tpq) = SSresults.SSnetdir_in + SSresults.SSnetdif_in + SSresults.SSnetLW_in; %[W/m^2/K] Net (All Ecocystem) Entropy in total energy fluxes reaching the ecosystem
  SSeco_net_out_store(1,tpq) = SSresults.SSnetdir_out + SSresults.SSnetdif_out + SSresults.SSnetLW_out + ...
                                    SSresults.SSLE_net + SSresults.SSH_net;%[W/m^2/K] Net (All Ecocystem) Entropy in total energy fluxes leaving the ecosystem
  SSeco_net_tot_store(1,tpq) = SSeco_net_in_store(1,tpq) - SSeco_net_out_store(1,tpq); %[W/m^2/K] Net (All Ecocystem) Entropy Budget 
  SSeco_net_ratio_store(1,tpq) = SSeco_net_out_store(1,tpq)/SSeco_net_in_store(1,tpq); %[W/m^2/K] Net (All Ecocystem) Entropy Ratio        
        
        % EFFECTIVE TEMPERATURES
  SSTl_net_store(1,tpq) = SSresults.Tl_net;                                 % [C] Effective Temperature (For All Ecosystem)  Based on Longwave
  SSTl_net2_store(1,tpq) = SSresults.Tl_net2;                               % [C] Effective Temperature (For All Ecosystem)  Based on MLCan model and averaging soil and canopy
  SSTeffent_store(1,tpq) = SSresults.Teffent;                               % [C] Effective Temperature (For All Ecosystem)  Based on entropy calculation
  SSXeffent_store(1,tpq) = SSresults.Xeffent;                               % [] Effective X_LW (For All Ecosystem) 

  % ENTROPY ASSOCIATED WITH WATER INFILTRATION  
  SH2O_store(:,tpq) =  SSresults.SH2O;                                      % [W/m2/K] Entropy due to water infiltration 
  EH2O_store(:,tpq) =  SSresults.EH2O;                                      % [W/m2] Energy associated with water budget in soil   
end

