function [ea,RH] = Ta_VPD_calculation(Ta,VPD)
% This function computes the atmospheric vapor pressure and also the
% relative humidity using atmospheric temperature and water vapor deficit

% Taken from An Introduction to Environmental Biophysics, p41
a = 0.611;  % [kPa]
%b = 17.27;
%c = 237.3;
b = 17.502;
c = 240.97; % [C]
esat = a*exp((b*Ta)./(c+Ta));
ea = esat - VPD; % [kPa]
RH = ea./esat; % [kPa]

