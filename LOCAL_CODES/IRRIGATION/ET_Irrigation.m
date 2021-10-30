function [IW, VARIABLES, CROP_GROWTH] = ET_Irrigation( VARIABLES, CROP_GROWTH, CONSTANTS, tt)

%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                              FUNCTION CODE                            %%
%%          Evapo-Transpiration based irrigation scheduling              %%


%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Rohit Nandan                                           %
%   Date         : February 12, 2020                                       %
%-------------------------------------------------------------------------%

%%

%BIOMASS_ESTIMATION Summary of this function goes here
%   Detailed explanation goes here


% ET_accu = VARIABLES.SOIL.ET_accu;
	
% For Irrigation Water ET: Rohit start
     
        tdt = 24 * 60 / CONSTANTS.timestep;

        tqt = tt + (CROP_GROWTH.DOY_start-1) * tdt;
        
       TR = VARIABLES.CANOPY.TR_can;
        
       E_soil = VARIABLES.SOIL.Esoil;
        
       ppt_gr = VARIABLES.SOIL.ppt_ground;
        
        plantingDate = CROP_GROWTH.plantingDate;
        harvestingDate = CROP_GROWTH.harvestingDate;

        
        Allow_Depl_Per = 0.50;
        
        Water_Holding_Capacity = 0.110;
        Rooting_Depth = 1.0 * 1000;
        
        Stored_Moisture = Water_Holding_Capacity * Rooting_Depth;
        
        Allow_Depl_Amount = Stored_Moisture * Allow_Depl_Per;
        
        Net_Irrigation = Allow_Depl_Amount;
        
        Total_ET = TR * 3600 + E_soil * 3600;
        
        
        if (tqt > plantingDate*tdt && tqt < harvestingDate*tdt)
           
            VARIABLES.SOIL.ET_accu = VARIABLES.SOIL.ET_accu + Total_ET - ppt_gr;
            
%             disp(VARIABLES.SOIL.ET_accu)
        
            if(VARIABLES.SOIL.ET_accu < 0)
               VARIABLES.SOIL.ET_accu = 0; 
            end
            
            if(VARIABLES.SOIL.ET_accu > Net_Irrigation)
        
%                 IW = Net_Irrigation/4;
                            
                CROP_GROWTH.time_delay=1;
                
                CROP_GROWTH.No_time=0;
                
%                 CROP_GROWTH.time_delay=CROP_GROWTH.time_delay+1;
                
%                 if (CROP_GROWTH.time_delay==3)
                    
                    VARIABLES.SOIL.ET_accu = VARIABLES.SOIL.ET_accu - Net_Irrigation;
                    
%                     CROP_GROWTH.time_delay=0;
                    
%                 end
        
            end
            
                 
            if (CROP_GROWTH.time_delay==1)
                
                IW = Net_Irrigation/14;
                
                CROP_GROWTH.No_time = CROP_GROWTH.No_time+1;
                
                if (CROP_GROWTH.No_time==14)
                CROP_GROWTH.time_delay=0;
                end
                
            else
            
                IW = 0;
                
            end
            
        
        else
            
            IW = 0;
            
        end
         
 
end
