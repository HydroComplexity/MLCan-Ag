function [ doy1, doy2, doy3, doy4, doy5, plantingDate, harvestingDate ] = Growing_Degree_Days(CROP_GROWTH, Max_Ta, Min_Ta, NoDay, TempCut, TempBase, daily_ppt, yy)

%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                              FUNCTION CODE                            %%
%%  Identified different growth stages using growing degree days and     %%
%%  precipitation based onset selection                                  %%

%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Rohit Nandan                                           %
%   Date         : November 18, 2019                                      %
%-------------------------------------------------------------------------%

%%

%GROWING_DEGREE_DAYS Summary of this function goes here
%   Detailed explanation goes here
% 
%     int NoDay = ForcingDaily->NoDay;
% 	double *AvgTa = new double[NoDay];
% 
% 	double *GDD = new double[Phenology->harvestingDate - Phenology->plantingDate];



% plantingDate = CROP_GROWTH.plantingDate;
harvestingDate = CROP_GROWTH.harvestingDate;

% GDD=zeros(harvestingDate - plantingDate);
GDD0=0;
GDD=0;

PPT_total=0;

plantingDate0 = 91;

doy5 = harvestingDate;

	for ii = 1 : NoDay
        
		if Max_Ta(ii) >= TempCut
			Max_Ta(ii) = TempCut;

        end
        
		if Min_Ta(ii) <= TempBase
			Min_Ta(ii) = TempBase;

        end
        
       AvgTa(ii) = (((Max_Ta(ii) + Min_Ta(ii)) / 2) - TempBase);
        
		if AvgTa(ii)<0
			AvgTa(ii) = 0;
        end
               
    end
    
    for ii = 1 : NoDay
        GDD0 = GDD0 + AvgTa (ii);
        
        if GDD0>90
                plantingDate0 = ii+ NoDay * (yy-1);
                break;
        end
    end
    
    for ii = plantingDate0 : 151
        PPT_total = PPT_total+ daily_ppt(ii);
            if PPT_total<50
                plantingDate0 = ii+ NoDay * (yy-1);
                break;
            end
    end
    
        
    if plantingDate0 < 91
        plantingDate = CROP_GROWTH.plantingDate;
    else
        plantingDate = plantingDate0;
    end
    
    
   
    disp ( ['Planting Date: ', num2str(plantingDate)]);
    
	for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
		GDD = GDD + AvgTa(ii);
		if GDD>CROP_GROWTH.GDDS1
			doy1 = ii + NoDay * (yy-1);
			disp ( ['Emmergence: ', num2str(doy1)]);
			GDD = 0;
			break;
        end
    end
        

	for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
		GDD = GDD + AvgTa (ii);
		if GDD>CROP_GROWTH.GDDS2
			doy2 = ii + NoDay * (yy-1);
            disp ( ['Iitial Vegetative stage: ', num2str(doy2)]);
			GDD = 0;
			break;
        end
    end
    
	for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
		GDD = GDD + AvgTa (ii);
		if GDD>CROP_GROWTH.GDDS3
			doy3 = ii + NoDay * (yy-1);
            disp ( ['Normal Vegetative stage: ', num2str(doy3)]);
			GDD = 0;
			break;
        end
    end

	for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
		GDD = GDD + AvgTa (ii);
		if GDD>CROP_GROWTH.GDDS4
			doy4 = ii + NoDay * (yy-1);
            disp ( ['Iitial Reproductive stage: ', num2str(doy4)]);
			GDD = 0;
			break;
        end
    end
    
	for ii = (plantingDate - NoDay * (yy-1)) : (harvestingDate - NoDay * (yy-1))
		GDD = GDD + AvgTa (ii);
		if GDD>CROP_GROWTH.GDDS5
			doy5 = ii + NoDay * (yy-1);
            disp ( ['Physical maturity stage: ', num2str(doy5)]);
			break;

        end
        
    end
    
    if harvestingDate > (doy5 + 20)
        harvestingDate = doy5 + 20;
    end
    
%     harvestingDate = CROP_GROWTH.harvestingDate;
    
    disp ( ['Harvesting Date: ', num2str(harvestingDate)]);
    
end

