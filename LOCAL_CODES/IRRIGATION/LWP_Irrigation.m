function [IW, VARIABLES, CROP_GROWTH] = LWP_Irrigation( VARIABLES, CROP_GROWTH, CONSTANTS, tt)

%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                              FUNCTION CODE                            %%
%%     Leaf water potential based irrigation scheduling                  %%


%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Rohit Nandan                                           %
%   Date         : February 12, 2020                                       %
%-------------------------------------------------------------------------%

%%

%LWP based irrigation Summary of this function goes here
%   Detailed explanation goes here


% For Irrigation Water Method 2: Rohit start

TR = VARIABLES.CANOPY.TR_can;

E_soil = VARIABLES.SOIL.Esoil;

ppt_gr = VARIABLES.SOIL.ppt_ground;

LWP = VARIABLES.CANOPY.psil_mean;

tdt = 24 * 60 / CONSTANTS.timestep;

tqt = tt + (CROP_GROWTH.DOY_start-1) * tdt;

plantingDate = CROP_GROWTH.plantingDate;
harvestingDate = CROP_GROWTH.harvestingDate;


Total_ET = TR * 3600 + E_soil * 3600;

if(rem(tqt, 24)== 0)
    CROP_GROWTH.New_T = tqt + (4:8);
end

if tqt > plantingDate*tdt && tqt < harvestingDate*tdt
    
    VARIABLES.SOIL.ET_accu = VARIABLES.SOIL.ET_accu + Total_ET - ppt_gr;
    
    if(VARIABLES.SOIL.ET_accu < 0)
        VARIABLES.SOIL.ET_accu = 0;
    end
    
    if(any(CROP_GROWTH.New_T(:)==tqt))
        
        % -0.6 is optimized
        if(LWP <= -0.6)
            
            CROP_GROWTH.time_delay=1;
            
            CROP_GROWTH.No_time=0;
            
            
            CROP_GROWTH.Net_Irrigation = VARIABLES.SOIL.ET_accu * 0.8;
            
            CROP_GROWTH.div = fix(CROP_GROWTH.Net_Irrigation/4);
            
            %                 CROP_GROWTH.div = fix(CROP_GROWTH.Net_Irrigation/12);
            
            %                     VARIABLES.SOIL.ET_accu = 0;
            %                 breakpoint=0;
            
            CROP_GROWTH.New_T=0;
            
            disp(['should be irrigated on ', num2str(tqt)])
            
        end
        
        %
    end
    
    
    if (CROP_GROWTH.time_delay==1)
        
        %                 Net_Irrigation = VARIABLES.SOIL.ET_accu * 0.8;
        %
        %                 div = fix(Net_Irrigation/12);
        
        if(CROP_GROWTH.div > 0)
            
            IW = CROP_GROWTH.Net_Irrigation/CROP_GROWTH.div;
            
            CROP_GROWTH.No_time = CROP_GROWTH.No_time + 1;
            
            if (CROP_GROWTH.No_time==CROP_GROWTH.div)
                CROP_GROWTH.time_delay=0;
                
                VARIABLES.SOIL.ET_accu = 0;
                
            end
            
        else
            
            IW=0;
            
        end
        
        
        
    else
        
        IW = 0;
        
    end
    
else
    
    IW = 0;
    
end


end
