function [ doy1, doy2, doy3, doy4, doy5, plantingDate, harvestingDate, CROP_GROWTH, VARIABLES ] =...
    Growing_Degree_Days(CROP_GROWTH, VARIABLES, CONSTANTS, PARAMS, SWITCHES, Max_Ta, Min_Ta, NoDay, TempCut, TempBase, yy, pp, tt)

%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                              FUNCTION CODE                            %%
%%  Identified different growth stages using growing degree days and Onset
%% using Soil moisture threshold

%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Rohit Nandan                                           %
%   Date         : November 18, 2019                                      %
%-------------------------------------------------------------------------%

%%


%GROWING_DEGREE_DAYS Summary of this function goes here
%   Detailed explanation goes here

doy1 = CROP_GROWTH.doy1(pp);
doy2 = CROP_GROWTH.doy2(pp);
doy3 = CROP_GROWTH.doy3(pp);
doy4 = CROP_GROWTH.doy4(pp);
doy5 = CROP_GROWTH.doy5(pp);

tdt = 24 * 60 / CONSTANTS.timestep;

tqt = tt + (CROP_GROWTH.DOY_start-1) * tdt;


plantingDate = CROP_GROWTH.plantingDate(pp);
harvestingDate = CROP_GROWTH.harvestingDate(pp);

ND=CROP_GROWTH.DOY_end-CROP_GROWTH.DOY_start;
ST=CROP_GROWTH.DOY_start;
ED=CROP_GROWTH.DOY_start;

% Vol_SM = mean(VARIABLES.SOIL.volliq(:));

if SWITCHES.Onset==1
    
    CROP_GROWTH.value=CROP_GROWTH.value+1;
    
    CROP_GROWTH.Vol_SM_all(CROP_GROWTH.value) = mean(VARIABLES.SOIL.volliq(:));
    
    if tqt>91*tdt && tqt<151*tdt
        
        if(rem(tqt, 24)== 0)
            
            CROP_GROWTH.Daily_SM = mean(CROP_GROWTH.Vol_SM_all);
            
            CROP_GROWTH.Vol_SM_all = null(24);
            
            CROP_GROWTH.value=0;
            
            
            %          if CROP_GROWTH.Daily_SM < 0.4
            if CROP_GROWTH.Daily_SM < CROP_GROWTH.Daily_SM_prev
                
                CROP_GROWTH.No_count=CROP_GROWTH.No_count+1;
                
                if CROP_GROWTH.No_count==10 && CROP_GROWTH.lock == 1
                    
                    plantingDate = (tqt/tdt) + NoDay * (yy-1);
                    
                    CROP_GROWTH.No_count=0;
                    
                    CROP_GROWTH.lock = 0;
                    
                end
                
            end
            %          end
            
            CROP_GROWTH.Daily_SM_prev=CROP_GROWTH.Daily_SM;
            
        end
        
    end
    
end

% CROP_GROWTH.plantingDate(pp) = plantingDate;

GDD=0;

for ii = 1 : ND
    
    if Max_Ta(ii) >= TempCut
        Max_Ta(ii) = TempCut;
    end
    if Min_Ta(ii) <= TempBase
        Min_Ta(ii) = TempBase;
    end
    
    AvgTa1(ii) = (((Max_Ta(ii) + Min_Ta(ii)) / 2) - TempBase);
    
    if AvgTa1(ii)<0
        AvgTa1(ii) = 0;
    end
    
end


BalanceB=zeros(1,ST-1);
BalanceE=zeros(1,NoDay-ED);

AvgTa=[BalanceB,AvgTa1,BalanceE];

if tqt == plantingDate * tdt
    
    doy5 = harvestingDate;
    
    disp ( ['Planting Date ', 'for Species ', num2str(pp),': ', num2str(plantingDate)]);
    
    for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
        GDD = GDD + AvgTa(ii);
        if GDD>PARAMS.CGM.GDDS1(pp)
            doy1 = ii + NoDay * (yy-1);
            disp ( ['Emmergence ', 'for Species ', num2str(pp),': ', num2str(doy1)]);
            GDD = 0;
            break;
        end
    end
    
    
    for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
        GDD = GDD + AvgTa (ii);
        if GDD>PARAMS.CGM.GDDS2(pp)
            doy2 = ii + NoDay * (yy-1);
            disp ( ['Iitial Vegetative stage ', 'for Species ', num2str(pp),': ', num2str(doy2)]);
            GDD = 0;
            break;
        end
    end
    
    for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
        GDD = GDD + AvgTa (ii);
        if GDD>PARAMS.CGM.GDDS3(pp)
            doy3 = ii + NoDay * (yy-1);
            disp ( ['Normal Vegetative stage ', 'for Species ', num2str(pp),': ', num2str(doy3)]);
            GDD = 0;
            break;
        end
    end
    
    for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
        GDD = GDD + AvgTa (ii);
        if GDD>PARAMS.CGM.GDDS4(pp)
            doy4 = ii + NoDay * (yy-1);
            disp ( ['Iitial Reproductive stage ', 'for Species ', num2str(pp),': ', num2str(doy4)]);
            GDD = 0;
            break;
        end
    end
    
    for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
        GDD = GDD + AvgTa (ii);
        if GDD>PARAMS.CGM.GDDS5(pp)
            doy5 = ii + NoDay * (yy-1);
            disp ( ['Physical maturity stage ', 'for Species ', num2str(pp),': ', num2str(doy5)]);
            break;
            
        end
        
    end
    
    if doy5>=harvestingDate || isnan(doy5) || doy5==0
        doy5=fix((harvestingDate + doy4)/2);
        disp ( ['Physical maturity stage ', 'for Species ', num2str(pp),': ', num2str(doy5)]);
    end
    
    if SWITCHES.Onset == 1
        
        if harvestingDate > (doy5 + 20)
            harvestingDate = doy5 + 20;
        end
        
    end
    
    
%     disp(num2str(doy5));
    
    disp ( ['Harvesting Date ', 'for Species ', num2str(pp),': ', num2str(harvestingDate)]);
    
%     
%     Lxi=(1: 24*NoDay)';
%     Lx=[1, 24*plantingDate,24*doy1,24*doy2,24*doy3,24*doy4,24*doy5,24*harvestingDate-1,24*harvestingDate, 24*NoDay]';
%     Ly=[0.1,0.1,0.1,3.5,5.6,3.5,2,1, 0.1,0.1]';
%     CROP_GROWTH.LAI_ref=interp1(Lx,Ly,Lxi,'pchip');
    
%     disp(CROP_GROWTH.LAI_ref(doy3*24));

    Lxi=(1: 24*NoDay)';
    Lx=[1, 24*plantingDate,24*doy1,24*doy2,24*doy3,24*doy4,24*doy5,24*harvestingDate-1,24*harvestingDate, 24*NoDay]';
    Ly=[0, 0 , 50, 600 , 850, 930, 1000, 920, 0, 0]';
    CROP_GROWTH.TotCmass_ref_s=interp1(Lx,Ly,Lxi,'pchip');
    
%     TCM_ref_h=0;
    CROP_GROWTH.TotCmass_ref=zeros(24*NoDay,1);
    for ii=2:24*NoDay
        CROP_GROWTH.TotCmass_ref(ii)=CROP_GROWTH.TotCmass_ref_s(ii)-CROP_GROWTH.TotCmass_ref_s(ii-1);
%         =TCM_ref_h;
    end
    
    CROP_GROWTH.TotCmass_ref(CROP_GROWTH.TotCmass_ref<0)=0;
    
end


