
function [ Max_Ta, Min_Ta, daily_ppt ] = Convert_to_Daily( Ta_tot, NoDay, PPT_tot, StepsT )

%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                              FUNCTION CODE                            %%
%% Convert the hourly temperature and precipitation into daily time      %%
%% interval forthe estimation of growing degree days                     %%

%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Rohit Nandan                                           %
%   Date         : November 18, 2019                                      %
%-------------------------------------------------------------------------%

%%

%CONVERT_TO_DAILY Summary of this function goes here
%   Detailed explanation goes here

td = 24 * 60 / StepsT;

for i = 1 : NoDay
    MaxDailyTa(i) = max(Ta_tot((1 + (i-1) * td) : (td + (i-1) * td)));
    MinDailyTa(i) = min(Ta_tot((1 + (i-1) * td) : (td + (i-1) * td)));
    Max_Ta(i) = MaxDailyTa(i);
    Min_Ta(i) = MinDailyTa(i);
    
    daily_ppt(i) = sum(PPT_tot((1 + (i-1) * td) : (td + (i-1) * td)));
end

% Ta_tot(:) = 0;

end

