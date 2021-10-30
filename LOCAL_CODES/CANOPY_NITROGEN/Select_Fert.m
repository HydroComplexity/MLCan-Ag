
for si=1:num_species
    
    VARIABLES.CANOPY.N_Demand(si) = (CROP_GROWTH.N_Demand(si)+0.032 * CROP_GROWTH.TotCmass_ref(tt-(NoDay_prev*24)))/2;
    
    VARIABLES.CANOPY.N_req(si) = VARIABLES.CANOPY.N_Demand(si) - (VARIABLES.UP_N_m2_sp(si))/24;
    
    file_name_Py   = sprintf('YrwdSpec_%sY%s.mat', num2str(si), num2str(Run_years(yy)));
    
    load('.\Temps\temporary.mat', 'config_path');
    
    fullpath_Py    =...
        fullfile(config_path, 'YearFr', file_name_Py);
    
    
    load(fullpath_Py, 'Fertilized', 'NitType', 'NitData')
    
    VARIABLES.Fertilized(si)=Fertilized;
    
    if Fertilized==1
        
        if (tt>=CROP_GROWTH.plantingDate(si)*24 && tt<=CROP_GROWTH.harvestingDate(si)*24)
            
            TotalSoilN=mean(VARIABLES.Amm(2:13, si).*roottr(:, si) +...
                VARIABLES.Nit(2:13, si).*roottr(:, si));
            
            VARIABLES.CANOPY.Net_Nreq(si)=(VARIABLES.CANOPY.Net_Nreq(si)+VARIABLES.CANOPY.N_req(si));
            
            VARIABLES.CANOPY.Net_Nreq(VARIABLES.CANOPY.Net_Nreq<=0)=0;
            
            
            if strcmpi(NitType, 'Manually')
                % Manually
                VARIABLES.Fert_Type(si)=1;
                
                Fert_DOY_list=cell2mat(NitData(:,1));
                %                 + NoDay_prev;
                %                 RunningDOY=tt/24;
%                 Famm=0; Fnit=0; Furea=0;
                
                if find(Fert_DOY_list==FORCING.doy)
                    
                    VARIABLES.CANOPY.N_Fert_DOY(si)=FORCING.doy;
                    
                    fr=find(Fert_DOY_list==FORCING.doy);
                    
                    Famm=cell2mat(NitData(fr,2));
                    Fnit=cell2mat(NitData(fr,3));
                    Furea=cell2mat(NitData(fr,4));
                    
                    
                    VARIABLES.N_Fert_amm(si)=Famm/24;
                    VARIABLES.N_Fert_nit(si)=Fnit/24;
                    VARIABLES.N_Fert_urea(si)=Furea/24;
                                        
                end
                
                
                VARIABLES.CANOPY.Net_Nreq(si)=0;
                
                VARIABLES.N_Fert_DOY_x(si)=VARIABLES.CANOPY.N_Fert_DOY(si);
                
                
            elseif strcmpi(NitType, 'Based on GS')
                % Crop growth
                VARIABLES.Fert_Type(si)=2;
                
                if (tt==CROP_GROWTH.plantingDate(si)*24)
                    
                    VARIABLES.CANOPY.N_Fert_DOY(si)=(tt/24)+1;
                    
                    Fert_req= 7;
                    
                    disp(['Fertilizer required here:', num2str(VARIABLES.CANOPY.N_Fert_DOY(si)), '|', num2str(Fert_req)])
                    
                elseif (tt==CROP_GROWTH.doy1(si)*24||tt==CROP_GROWTH.doy2(si)*24||...
                        tt==CROP_GROWTH.doy3(si)*24||tt==CROP_GROWTH.doy4(si)*24)
                    
                    VARIABLES.CANOPY.N_Fert_DOY(si)=(tt/24)+1;
                    
                    Fert_req=VARIABLES.CANOPY.Net_Nreq(si)-TotalSoilN;
                    
                    disp(['Fertilizer required here:', num2str(VARIABLES.CANOPY.N_Fert_DOY), '|', num2str(Fert_req)])
                    
                    Fert_req(Fert_req<=0)=0;
                    
                else
                    
                    VARIABLES.CANOPY.N_Fert_DOY(si)=0;
                    
                    Fert_req=0;
                    
                end
                
            elseif strcmpi(NitType, 'With Irrigation')
                % Based on Irrigation
                VARIABLES.Fert_Type(si)=3;
                
                if (tt==CROP_GROWTH.plantingDate(si)*24)
                    VARIABLES.CANOPY.N_Fert_DOY(si)=(tt/24)+1;
                    Fert_req= 7;
                    disp(['Fertilizer required here:', num2str(VARIABLES.CANOPY.N_Fert_DOY(si)), '|', num2str(Fert_req)])
                elseif(IWD==1)
                    VARIABLES.CANOPY.N_Fert_DOY(si)=(tt/24)+1;
                    Fert_req=VARIABLES.CANOPY.Net_Nreq(si)-TotalSoilN;
                    disp(['Fertilizer required here:', num2str(VARIABLES.CANOPY.N_Fert_DOY(si)), '|', num2str(Fert_req)])
                    Fert_req(Fert_req<=0)=0;
                else
                    VARIABLES.CANOPY.N_Fert_DOY(si)=0;
                    Fert_req=0;
                end
                
                
            elseif strcmpi(NitType, 'Leaf Nitrogen')
                
                VARIABLES.Fert_Type(si)=4;
                
                Nleaf_ga=sum(VARIABLES.NCan.Nleaf(:,si))*CROP_GROWTH.LAIsim_can(si);
                
                Nleaf_thresh=Nleaf_ga-VARIABLES.NCan.Nleaf_ga(si);
                
                if (Nleaf_thresh<0.5)
                    VARIABLES.CANOPY.N_Fert_DOY(si)=(tt/24)+1;
                    Fert_req=VARIABLES.CANOPY.Net_Nreq(si)-TotalSoilN;
                    disp(['Fertilizer required here:', num2str(VARIABLES.CANOPY.N_Fert_DOY(si)), '|', num2str(Fert_req)])
                    Fert_req(Fert_req<=0)=0;
                end
                
                VARIABLES.NCan.Nleaf_ga(si)=Nleaf_ga;
                
            end
            
        else
            
            VARIABLES.CANOPY.N_Fert_DOY(si)=0;
            Fert_req=0;
            
        end
        
    else
        VARIABLES.CANOPY.N_Fert_DOY(si)=0;
        
        Fert_req=0;
    end
    
    
    if VARIABLES.Fert_Type(si)>1
        
        VARIABLES.CANOPY.Net_Nreq(si)=0;
        
        VARIABLES.N_Fert_DOY_x(si)=VARIABLES.CANOPY.N_Fert_DOY(si);
        
        VARIABLES.N_Fert_amm(si)=0.25*Fert_req/24;
        VARIABLES.N_Fert_nit(si)=0.25*Fert_req/24;
        VARIABLES.N_Fert_urea(si)=0.50*Fert_req/24;
        
    end
    
end
