function [zen] = ZEN_calculation(LAT,LONG,hours,doys)
%
%   Compute the zenith angle for a particular site and range of times
%   
%   Taken From Campbell and Norman, Eqns 11.1 - 11.5
%   (An Introduction to Environmental Biophysics, p168-170)
%   Written By : Darren Drewry
%   Modified By: Dongkook Woo
%
%   INPUTS:
%       LAT - latitude of site [degrees]
%       hours - decimal hours (0..24)
%       doys - integer day of year
%       LONG   - Longitude [degrees]
%
%  OUTPUTS:
%       zen  - Zenith angle [degrees]
%       dec  - Declination angle [degrees]
%

d2r = pi/180;   % Conversion factor from degrees to radians
r2d = 180/pi;   % Conversion factor from radians to degrees

% Dongkook Woo - Edit
% doy should start with 1
doys=floor(round((doys+1).*10^10)./10^10);
% Dongkook Woo - Edit End

% Solar Declination Angle [rad] - ranges from -23.45 degrees to + 23.45 degrees
d1 = d2r*(356.6 + 0.9856*doys); % [rad]
d2 = sin(d1);                   % [rad]
d3 = 1.9165 * d2;               % [deg]
d4 = 278.97 + 0.9856*doys + d3; % [deg]
d5 = sin( d2r*d4 );             % [rad]
d6 = 0.39785 * d5;              % [deg]
dec = asin( d6 );

ff = d2r*(279.575 + 0.9856*doys);    % [rad]
ET = (-104.7*sin(ff) + 596.2*sin(2*ff) + 4.3*sin(3*ff) - 12.7*sin(4*ff) - 429.3*cos(ff) - 2.0*cos(2*ff) + 19.3*cos(3*ff)) / 3600;

% Compute the longitude correction 1/15 grades
% Dongkook Woo - Edit
%LC=LONG*(1/15);
n=ceil(LONG./15);
LC=(15*n-LONG)*(1/15);
% Dongkook Woo - Edit End

t0 = 12 - LC - ET;              %[hours] 

LAT = d2r*LAT;  % [rad]
zen = acos( sin(LAT)*sin(dec) + cos(LAT)*cos(dec).*cos(d2r*15*(hours-t0)) );    %[Zenith Angle]

zen = zen * r2d;  % [degrees]
dec = dec * r2d;  % [degrees]