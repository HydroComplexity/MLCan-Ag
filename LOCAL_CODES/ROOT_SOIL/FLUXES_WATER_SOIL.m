
function [VARIABLES] = FLUXES_WATER_SOIL (PARAMS, VARIABLES, CONSTANTS,...
    FORCING, SWITCHES, si)

%=========================================================================
% Solve surface energy balance (With snow and No Litter)
%
% Written by Juan Quijano, UIUC, 2013
% All rights reserved!
%
%------------------------- Input Variables -------------------------------
%       PARAMS          % PARAMS structure
%       VARIABLES       % VARIABLES structure
%       CONSTANTS       % CONSTANTS structure
%       FORCING         % CONSTANTS structure
%       SWITCHES        % SWITCHES structure
%------------------------- Output Variables ------------------------------
%       VARIABLES       % VARIABLES structure
% 
%========================================================================
    
    % De reference block
    %*********************************************************************
    % PARAMS
    thetals = PARAMS.Soil.thetals(si);                                          % Litter Soil Moisture at saturation
    thetafc = PARAMS.Soil.thetafc(si);                                          % Litter moisture at field capacity
    km = PARAMS.Soil.km(si);                                                    % parameter to compute drainage from litter
    bm = PARAMS.Soil.bm(si);                                                    % parameter to compute drainage from litter
    Tcr_atm = PARAMS.Tcr_atm(si);                                               % [K]Critical Atmospheric Temperature around 273 K
    Tf = PARAMS.Soil.Tf(si);                                                    % [K] Freezing Temperature of Fresh Water 
    rho_liq = PARAMS.Soil.rho_liq(si);                                          % [kg / m^3] Density of Liquid Water
    rho_ice = PARAMS.Soil.rho_ice(si);                                          % [km / m^2] Density of ice
    HC_snow = PARAMS.Soil.HC_snow(si);                                          % [] Maximum fraction water capacity that snow can hold 
    rhosn_min = PARAMS.Soil.rhosn_min(si);                                      % [kg / m^3] Minimum Snow Density     
    fcomp = PARAMS.Soil.fcomp(si);                                              % Compactation parameter          
        
    % VARIABLES
    wliqsl = VARIABLES.SOIL.wliqsl(si);                                         % [kg / m^2] Liquid Water Density per unit area 
    wicesl = VARIABLES.SOIL.wicesl(si);                                         % [kg / m^2] Ice Water Density per unit area
    wsn = VARIABLES.SOIL.wsn(si);                                               % [kg / m^2] Total Snow Density per unit area
    zliqsl = VARIABLES.SOIL.zliqsl(si);                                         % [mm] Liquid water depth
    zicesl = VARIABLES.SOIL.zicesl(si);                                         % [mm] Ice water depth
    zsn = VARIABLES.SOIL.zsn(si);                                               % [mm] Snow water depth
    deltawice = VARIABLES.SOIL.deltawice(si);                                   % [] Change in volumetric volume of ice
    rhosn = VARIABLES.SOIL.rhosn(si);                                           % [kg / m^3] Density of snow
    Tli = VARIABLES.SOIL.Tli(si);                                               % [C]Temperature of snow pack    
    ppt_ground = VARIABLES.SOIL.ppt_ground(si);                                 % [mm] precipitation reaching ground
    snow_tcount = VARIABLES.SOIL.snow_tcount(si);                               % [] Number of time steps to account snow age    
    dzlit_mm = VARIABLES.SOIL.dzlit_m(si)*1000;                                 % Thickness of litter [mm]
    LEsl = VARIABLES.SOIL.LEl(si);                                              % [W/m^2] Latent Heat from snow pack  
    LEs = VARIABLES.SOIL.LEs(si);                                               % [W/m^2] Latent Heat from the soil 

    % CONSTANTS
    grav = CONSTANTS.grav;                                                  % [m / s^2] Gravity Acceleration   
    dtime = CONSTANTS.dtime;                                                % [s] Time Step
    Lv_kg = CONSTANTS.Lv_kg;                                                % [J/kg] Latent Heat of Evaporation 
    Lf_kg = CONSTANTS.Lv_kg;                                                % [J/kg] Latent Heat of Fusion
    Ls_kg = Lv_kg + Lf_kg;                                                  % [J/kg] Latent Heat of Sublimation 
    
    % FORCING
    Ta = FORCING.Ta;                                                        % [C] Atmospheric temperature
    
    %*********************************************************************
    %i. Compute LE fluxes 
    Esoil = LEs / Lv_kg;                                                    % Evaporation From Soil [mm/s]
    if zicesl > 0;         
        Esl_ice = LEsl / Ls_kg * (1000/rho_ice);                           %[mm ice / s]Sublimation From Sow Pack 
        Esl = Esl_ice;
        Esl_liq = 0;
    else        %(Very likely Esl = 0)
        Esl_ice = 0;                  
        Esl_liq = LEsl / (Lv_kg) * (1000/rho_liq);                         %[mm liq / s] Evaporation From Sow Pack 
        Esl = Esl_liq;        
    end 
        
    % ii. snow density of snow reaching ground [kg/m3]
    if Ta < -15 
        rhosn_new = 50;
    elseif (Ta >= -15 && Ta <= 2)
        rhosn_new = 50 + 1.7*(Ta+15)^(1.5); 
    elseif Ta > 2
        rhosn_new = 50 + 1.7*(17)^(1.5);
    end    
    % iii. compute rainfall
    pptrate_ground = ppt_ground / dtime; 
    if Ta < (Tcr_atm - Tf)
        zice_gro = pptrate_ground*dtime *rho_liq/rho_ice;                  %[mm ice] 
        zliq_gro = 0;
        snow_tcount = 0;
    else        
        zice_gro = 0;
        zliq_gro = pptrate_ground*dtime;
        rhosn_new = rho_liq;
        snow_tcount = snow_tcount + 1;
    end
        
    % iv. Computation of Drainage
    if (zicesl == 0 && Ta > (Tcr_atm - Tf))                                % No snow, rain water infiltrates in the soil directly
            drainlitter = 0;
            zliq_gro = 0;
            if SWITCHES.soilevap   
               qinfl = pptrate_ground - Esoil;                             %[mm/s] Infiltration into soil with extraction of soil evaporation [mm]
            else
               qinfl = pptrate_ground;                                     %[mm/s]  
            end
            qinflL = 0;
            net_qinflL = 0;
    else                   
            drainlitter = max((wliqsl/wsn - HC_snow)*wsn/(rho_liq)*1000,0); % Water that drains into soil from litter [mm]
            % v. Compute infiltration into snow-litter pack and into soil    
            if SWITCHES.soilevap   
                qinfl = drainlitter / dtime - Esoil;                       %[mm/s] Infiltration into soil with extraction of soil evaporation [mm]
            else
                qinfl = drainlitter / dtime;                               %[mm/s]  
            end
            qinflL = pptrate_ground - Esl;                                 %[mm/s] Water that infiltrates into litter    
            net_qinflL = pptrate_ground - Esl - drainlitter/dtime;         %[mm/s] Net Water that infiltrates into litter (with no evaporation no drainage)
    end    
    % vi. Compute mass balance
    %remains 
    remainsl = wliqsl + wicesl - (drainlitter/1000).*rho_liq ...
                - (Esl_liq/1000)*dtime*rho_liq - (Esl_ice/1000)*dtime*rho_ice; %[kg]
    % liq and ice
    dwliq_mm = -deltawice/rho_liq*1000;                                             %[mm liq]
    zliqsl_new = zliqsl + (zliq_gro - drainlitter - Esl_liq*dtime + dwliq_mm);      %[mm liq]    
    if zliqsl_new < 0 
        qinfl = qinfl - (0 - zliqsl_new)/dtime;
        drainlitter = min(drainlitter - (0 - zliqsl_new),0);                        %[mm liq]
        zliqsl_new = 0;                                                             %[mm liq]
        % correct qinfl
    end
        
    wliqsl_new = (zliqsl_new/1000)*rho_liq;                                         %[kg H2O /m2]
        
    dwice_mm = deltawice/rho_ice*1000;                                              %[mm ice]
    zicesl_new = zicesl + (zice_gro - Esl_ice*dtime + dwice_mm);                    %[mm ice]    
    if (zicesl_new < 0 && SWITCHES.soilevap)
        qinfl = qinfl + zicesl_new/dtime;                                           %[mm ice]
    end    
    zicesl_new = max(0,zicesl_new);                                                 %[mm ice]
    wicesl_new = (zicesl_new/1000)*rho_ice;                                         %[kg H2O /m2] 
         
    
    %vii. Compute snow density
    % ***** compute CRm ***** 
    if rhosn <= 150 
        c3 = 1;
    else
        c3 = exp(-0.046*(rhosn-150));
    end
    if zliqsl == 0
        c4 = 1;
    else
        c4 = 2;
    end
    CRm = 2.778*10^(-6)*c3*c4*exp(-0.04*(-Tli)); 
    % ***** compute CRo ***** 
    c5 = 0.08;
    c6 = 0.021;
    Ps = 0.5*grav*rho_liq*(zice_gro + zliq_gro + fcomp*zsn);
    CRo = (Ps/(3.6*10^6))*exp(-c5*(-Tli))*exp(-c6*rhosn);
    % ***** compute change *****
    drhosn = (CRm + CRo)*rhosn*dtime; 
    rhosn_updt = rhosn + drhosn;
    if (zicesl < 0.01*dzlit_mm*thetals || rhosn_updt > rho_liq)                         
        rhosn_updt = rho_liq;                                              % rho_updt is equal to rho_liq if ice is very small
    end
    if rhosn_updt < rhosn_min
        rhosn_updt = rhosn_min;
    end    
    %compact_net = rhosn_updt/rhosn;
    
    % ponderate rho of water in snow-litter pack.    
    % The final snow density is computed by a mass wighted average of the
    % different densities involved
    rhosn = ((rhosn_updt*wsn) - ((Esl_liq/1000)*dtime*rho_liq)*rho_liq ... 
               - ((Esl_ice/1000)*dtime*rho_ice)*rho_ice - ((drainlitter/1000).*rho_liq)*rho_liq ...
               + ((zice_gro/1000)*rho_liq)*rhosn_new + ((zliq_gro/1000)*rho_liq)*rhosn_new) / ...
                 (remainsl + (zice_gro/1000)*rho_liq + (zliq_gro/1000)*rho_liq);
    rhosn = max(50,rhosn);                                                 % Snow density can not be lower than 50

    if zicesl < 0.01*dzlit_mm*thetals                         
        rhosn = rho_liq;                                                   % rho is equal to rho_liq if ice is very small
    end              
               
    % vii. Save new variables.
    % rhosn
    VARIABLES.SOIL.rhosn(si) = rhosn;                                          % [kg / m^3] Density of snow
    % w
    VARIABLES.SOIL.wliqsl(si) = wliqsl_new;                                    % [kg / m^2] Liquid Water Density per unit area  
    VARIABLES.SOIL.wicesl(si) = wicesl_new;                                    % [kg / m^2] Ice Water Density per unit area  
    % z
    VARIABLES.SOIL.zicesl(si) = zicesl_new;                                    % [mm] Ice water depth 
    VARIABLES.SOIL.zliqsl(si) = zliqsl_new;                                    % [mm] Liquid water depth  
    VARIABLES.SOIL.zicesl_prev(si) = zicesl;                                    % [mm] Ice water depth in previous time step
    VARIABLES.SOIL.zliqsl_prev(si) = zliqsl;                                    % [mm] Liquid water depth in previous time step

    if wicesl_new > 0     
        wsn = wicesl_new + wliqsl_new;                                      % [kg / m^2]
        zsn = (wsn/rhosn)*1000;                                             % [mm] Snow depth
    else
        zsn = 0;
        wsn = 0;
    end
    
    VARIABLES.SOIL.wsn(si) = wsn;                                               % [kg / m^2] Snow Water Density per unit area     
    VARIABLES.SOIL.zsn(si) = zsn;                                               % [mm]
        
    % vols
    VARIABLES.SOIL.volliqli(si) = nan;                                         % [] Volumetric Water Litter Content     
    VARIABLES.SOIL.voliceli(si) = nan;                                         % [] Volumtetric Ice Litter Content    
    VARIABLES.SOIL.voltotli(si) = nan;                                         % [] Volumetric Ice and Water Litter Content    

    VARIABLES.SOIL.volliqsn(si) = zliqsl_new/zsn;                              % [] Volumetric Liquid Water in Snow Content     
    VARIABLES.SOIL.volicesn(si) = zicesl_new/zsn;                              % [] Volumetric Ice Water in Snow Conent       
    VARIABLES.SOIL.voltotsn(si) = (zliqsl_new + zicesl_new)/zsn;               % [] Volumetric Ice and Water Snow Content     
    
    %Other variables 
    VARIABLES.SOIL.CRm(si) = CRm;                                              % Snow Compaction Factor 
    VARIABLES.SOIL.CRo(si) = CRo;                                              % Snow Compaction Factor 
    VARIABLES.SOIL.drhosn(si) = drhosn;                                        % [kg / m^3] Change in snow density  
    
    VARIABLES.SOIL.qinflL(si) = qinflL;                                        %[mm/s] Water that infiltrates into snow  
    VARIABLES.SOIL.net_qinflL(si) = net_qinflL;                                %[mm/s] Net water that infiltrates into snow excluding drainage in soil    
    VARIABLES.SOIL.qinfl(si) = qinfl;                                          %[mm/s] Net Infiltration into the soil  
    VARIABLES.SOIL.drainlitter(si) = drainlitter;                              %[mm] Drainage from snow into soil  
    VARIABLES.SOIL.pptrate_ground(si) =  pptrate_ground;                       %[mm/s] Precipitation reaching the ground 
    VARIABLES.SOIL.Esoil(si) = Esoil;                                          %[mm/s] Evaporation soil  
    VARIABLES.SOIL.Esl(si) = Esl;                                              %[mm/s] Sublimation from snow 
    VARIABLES.SOIL.zliq_gro(si) = zliq_gro;                                    %[mm] Depth of liquid water reaching snow 
    VARIABLES.SOIL.zice_gro(si) = zice_gro;                                    %[mm] Depth of ice water reaching snow  
    %
    VARIABLES.SOIL.Esl_ice(si) = Esl;                                          %[mm/s] Sublimation from ice in the snow  
    VARIABLES.SOIL.Esl_liq(si) = nan;                                          %[mm/s] Evaporation from snow (assumed zero)  
    % Count time since last snow
    VARIABLES.SOIL.snow_tcount(si) = snow_tcount;                              %[] Number of time steps with snow 
