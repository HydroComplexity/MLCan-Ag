
for si=1:num_species
    
    Tl_canopy_prof(:,si) = (nansum(Tl_sun(:,si).*fLAIz(:,si),2)).*fsun(:) + (nansum(Tl_shade(:,si).*fLAIz(:,si),2)).*fshade(:);
    Tl_mean_all = sum(Tl_canopy_prof(vinds,si).*LADnorm_all(vinds,si))/sum(LADnorm_all(:,si));
    
    Tl_canopy_prof_up(:,si) = (nansum(Tl_sun_up(:,si).*fLAIz(:,si),2)).*fsun(:) + (nansum(Tl_shade_up(:,si).*fLAIz(:,si),2)).*fshade(:);
    Tl_mean_up = sum(Tl_canopy_prof_up(vinds,si).*LADnorm_all(vinds,si))/sum(LADnorm_all(:,si));
    
    Tl_canopy_prof_lo(:,si) = (nansum(Tl_sun_lo(:,si).*fLAIz(:,si),2)).*fsun(:) + (nansum(Tl_shade_lo(:,si).*fLAIz(:,si),2)).*fshade(:);
    Tl_mean_lo = sum(Tl_canopy_prof_lo(vinds,si).*LADnorm_all(vinds,si))/sum(LADnorm_all(:,si));
    
    if isnan(Tl_mean_up)
        Tl_mean_up=VARIABLES.CANOPY.TAz(1,si);
    end
    
    if isnan(Tl_mean_lo)
        Tl_mean_lo=VARIABLES.CANOPY.TAz(1,si);
    end
    
    if isnan(Tl_mean_all)
        Tl_mean_all=VARIABLES.CANOPY.TAz(1,si);
    end
    
    VARIABLES.CANOPY.Tl_mean_up(si) = Tl_mean_up;
    VARIABLES.CANOPY.Tl_mean_lo(si) = Tl_mean_lo;
    VARIABLES.CANOPY.Tl_mean(si) = Tl_mean_all;
    
    
    if (SWITCHES.Irr == 1)
        [IW, VARIABLES, CROP_GROWTH] = ET_Irrigation(VARIABLES, CROP_GROWTH, CONSTANTS, tt);
        CROP_GROWTH.CWSI(si) = 0;
        
    elseif (SWITCHES.Irr == 2)
        [IW, VARIABLES, CROP_GROWTH] = SC_Irrigation(VARIABLES, CROP_GROWTH, CONSTANTS, tt);
        CROP_GROWTH.CWSI(si) = 0;
        
    elseif (SWITCHES.Irr == 3)
        [IW, VARIABLES, CROP_GROWTH] = CWSI_Irrigation(VARIABLES, CROP_GROWTH, CONSTANTS, tt);
        
    elseif (SWITCHES.Irr == 4)
        [IW, VARIABLES, CROP_GROWTH] = LWP_Irrigation(VARIABLES, CROP_GROWTH, CONSTANTS, tt);
        CROP_GROWTH.CWSI(si) = 0;
        
    else
        IW = 0;
        CROP_GROWTH.CWSI(si) = 0;
        
    end
    
    IWD(si)=0;
    
    if(IW==0 && VARIABLES.SOIL.IW(si)>0)
        IWD(si)=1;
    end
    
    VARIABLES.SOIL.IW(si) = IW;
    
end