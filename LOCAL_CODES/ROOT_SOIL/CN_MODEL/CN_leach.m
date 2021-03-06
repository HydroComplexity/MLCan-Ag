function [LCH_CC, TLCH_CC] = CN_leach(PARAMS, SWITCHES, CC, aa, qq, sm)

nspecies = PARAMS.CanStruc.nspecies;

% CC Concentration variable [gr/m3]
% aa solubility of element in study
if SWITCHES.CN.Bioturbation
    nl_soil=PARAMS.nl_soil+1;%       nl_soil = # soil layers
else
    nl_soil=PARAMS.nl_soil;%       nl_soil = # soil layers
end

leach_on = SWITCHES.CN.leach_on;%       leach_on = Leaching ON or OF [1 0]

% Declare matrix
LCH_CC = zeros (nl_soil,nspecies);

for si=1:nspecies
    for ii = 1:nl_soil
        %   positive value = loss (subtracted from budget)
        %   NOTE: first H2O flux value is from infiltration into soil column
        LCH_CC(ii,si) = 0;
        if (leach_on)
            if (ii>1)   % transport through top layer interface
                if (qq(ii)>0)  % downward transport into cell
                    LCH_CC(ii,si) = LCH_CC(ii,si) - (aa/sm((ii-1),si)*CC((ii-1),si)*qq(ii,si));% [gr/m^2/d]
                    if isinf(LCH_CC(ii,si))
                        LCH_CC(ii,si)=0;
                    end
                else           % upward transport out of cell
                    LCH_CC(ii,si) = LCH_CC(ii,si) - (aa/sm(ii,si)*CC(ii,si)*qq(ii,si));% [gr/m^2/d]
                    if isinf(LCH_CC(ii,si))
                        LCH_CC(ii,si)=0;
                    end
                end
            end
            % transport through bottom layer interface
            if (qq((ii+1),si)>0)  % downward transport out of cell
                LCH_CC(ii,si) = LCH_CC(ii,si) + (aa/sm(ii,si)*CC(ii,si)*qq((ii+1),si));% [gr/m^2/d]
                if isinf(LCH_CC(ii,si))
                    LCH_CC(ii)=0;
                end
            else           % upward transport into cell
                LCH_CC(ii,si) = LCH_CC(ii,si) + (aa/sm((ii+1),si)*CC((ii+1),si)*qq((ii+1),si));% [gr/m^2/d]
                if isinf(LCH_CC(ii,si))
                    LCH_CC(ii,si)=0;
                end
            end
            if (ii==nl_soil)
                TLCH_CC(si) = (aa/sm(ii,si)*CC(ii,si)*qq((ii+1),si));% [gr/m^2/d]
            end
        end
    end
end

%     TLCH_CC = TLCH_CC(:);
%     LCH_CC = LCH_CC(:);

CCnew = CC - LCH_CC;
ind = CCnew < 0;
CCnew(ind) = 0;

%     PPPPPPPP;

