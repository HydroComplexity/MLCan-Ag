function [dwat, smp, hk, smp_weight, thsatfrac_weight, q_save] = ...
    SOILMOISTURE(SWITCHES, VERTSTRUC, PARAMS, VARIABLES, CONSTANTS, si)

%=========================================================================
% This code solves the model for water flow in the soil system, i.e.,
% the Richards equation with a sink term. The discretization results in
% a tridiagonal system of equation. The upper boundary condition is set
% to the infiltration from rain, and the lower boundary is set to the
% hydraulic conductivity of the bottom soil layer.
%
%------------------------- Input Variables -------------------------------
%         nl_soil       % number of soil layers
%         dtime         % time step [s]
%         wimp          % water impremeable if porosity less than wimp
%         smpmin        % restriction for min of soil poten. [mm]
%         psi0          % saturated soil suction [mm]
%         bsw           % Clapp-Hornberger "B"
%         porsl         % saturated volumetric soil water content(porosity)
%         hksati        % hydraulic conductivity at saturation [mm/s]
%         z             % layer depth [mm]
%         dz            % layer thickness [mm]
%         zi            % interface level below a "z" level [mm]
%         qinfl         % net water input from top [mm/s]
%         TR           % actual transpiration [mm/s]
%         rootfr        % root resistance of a layer, all layers add to 1.0
%         eff_porosity  % effective porosity of soil
%         vol_liq       % soil water per unit volume [mm/mm]
%         volice        % soil water per unit volume [mm/mm]
%         hr            % option for hydraulic redistribution
%
%------------------------- Output Variables ------------------------------
%         dwat          % change of soil water [m3/m3]
%         smp           % soil matrix potential [mm]
%         hk            % soil hydraulic conductivity [mm/s]
%         smp_weight    % smp weighted by root uptake [MPa]
%
%-------------------------- local variables ------------------------------
%        amx            % "a" left off diagonal of tridiagonal matrix
%        bmx            % "b" diagonal column for tridiagonal matrix
%        cmx            % "c" right off diagonal tridiagonal matrix
%        rmx            % "r" forcing term of tridiagonal matrix
%
%=========================================================================
%
% Determine the soil hydraulic conductivity and its derivative
%-------------------------------------------------------------------------
% Set the hydraulic conductivity to zero if effective porosity is < 5%
% in any of two neighbour layers or if liquid water content (theta) < 0.001


%*************************************************************************
%                          DE-REFERENCE BLOCK
%*************************************************************************

rhc = SWITCHES.rhc;
hr = SWITCHES.HR_on;

TR = VARIABLES.CANOPY.TR_can;

smp1 = VARIABLES.SOIL.smp(:,si);
volliq = VARIABLES.SOIL.volliq(:,si);
qinfl = VARIABLES.SOIL.qinfl(si);

rpp = VARIABLES.ROOT.rpp(:,si);
krad = VARIABLES.ROOT.krad(:,si);

z = VERTSTRUC.znsmm;
dz = VERTSTRUC.dzsmm;
zi = VERTSTRUC.zhsmm;
rootfr = VERTSTRUC.rootfr(:,si);
thetadry = VERTSTRUC.theta_dry;
eff_porosity = VERTSTRUC.eff_poros;
bsw = VERTSTRUC.bsw;
hksati = VERTSTRUC.HKsat;
porsl = VERTSTRUC.porsl;
psi0 = VERTSTRUC.psi0;

nl_soil = PARAMS.Soil.nl_soil;
smpmin = PARAMS.Soil.smpmin;
wimp = PARAMS.Soil.wimp;

dtime = CONSTANTS.dtime;
mmH2OtoMPa = CONSTANTS.mmH2OtoMPa;
%*************************************************************************
%*************************************************************************

for j = 1:nl_soil
    if((eff_porosity(j) < wimp) ...
            | (eff_porosity(min(nl_soil,j+1)) < wimp) ...
            | (volliq(j) <= 1.e-3))
        hk(j) = 0;
        dhkdw(j) = 0;
    else
        s1 = 0.5*(volliq(j)/porsl(j)+volliq(min(nl_soil,j+1))...
            /porsl(min(nl_soil,j+1)));
        s2 = hksati(j)*s1^(2*bsw(j)+2);
        hk(j) = s1*s2;
        dhkdw(j) = (2*bsw(j)+3)*s2*0.5/porsl(j);
        if(j == nl_soil)
            dhkdw(j) = dhkdw(j) * 2;
        end
    end
end

% Determine the soil matrix potential and its derivative
%-------------------------------------------------------------------------
for j = 1:nl_soil
    if(porsl(j) < 1.e-6)     % bed rock
        s_node = 0.001;
        smp1(j) = psi0(j);
        dsmpdw(j) = 0;
    else
        s_node = min(1,volliq(j)/porsl(j));
        s_node = max(s_node, 0.001);
        smp1(j) = psi0(j)*s_node^(-bsw(j));
        smp1(j) = max(smpmin, smp1(j));        % limit soil suction
        dsmpdw(j) = -bsw(j)*smp1(j)/(s_node*porsl(j));
    end
end

% Setup the vectors, 'a', 'b', 'c', and 'r' for tridiagonal matrix solution
%-------------------------------------------------------------------------

if hr == 1   % if hydraulic redistribution is allowed to take place
    
    % For the top soil layer
    j      = 1;
    qin    = qinfl;
    
    den    = (z(j+1)-z(j));
    num    = (smp1(j+1)-smp1(j)) - den;
    qout   = -hk(j)*num/den;
    dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den;
    dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den;
    
    rmx(j) =  qin - qout - krad(j)*(smp1(j)-rpp(j));
    amx(j) =  0;
    bmx(j) =  dz(j)/dtime + dqodw1;
    cmx(j) =  dqodw2;
    
    q_save(1) = qin;
    
    % For the middile soil layers
    for j = 2:nl_soil - 1
        den    = (z(j) - z(j-1));
        num    = (smp1(j)-smp1(j-1)) - den;
        qin    = -hk(j-1)*num/den;
        dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den;
        dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den;
        
        den    = (z(j+1)-z(j));
        num    = (smp1(j+1)-smp1(j)) - den;
        qout   = -hk(j)*num/den;
        dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den;
        dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den;
        
        rmx(j) =  qin - qout - krad(j)*(smp1(j)-rpp(j));
        amx(j) = -dqidw0;
        bmx(j) =  dz(j)/dtime - dqidw1 + dqodw1;
        cmx(j) =  dqodw2;
        
        q_save(j) = qin;
    end
    
    % For the bottom soil layer
    j      = nl_soil;
    den    = (z(j) - z(j-1));
    num    = (smp1(j)-smp1(j-1)) - den;
    qin    = -hk(j-1)*num/den;
    dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den;
    dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den;
    
    qout   =  hk(j);
    dqodw1 =  dhkdw(j);
    
    rmx(j) =  qin - qout - krad(j)*(smp1(j)-rpp(j));
    amx(j) = -dqidw0;
    bmx(j) =  dz(j)/dtime - dqidw1 + dqodw1;
    cmx(j) =  0;
    
    q_save(nl_soil) = qin;
    q_save(nl_soil+1) = qout;
    
    
    smp1 = smp1(:);
    rpp = rpp(:);
    
    % Determine levels with and without HR
    noHRinds = find(smp1>rpp);
    HRinds = find(smp1<=rpp);
    
    % H2O Uptake from each level
    UPz(noHRinds) = krad(noHRinds).*(smp1(noHRinds)-rpp(noHRinds));
    UPz(HRinds) = 0;
    
    % H2O Redistributed into each level
    HRz(HRinds) = krad(HRinds).*(smp1(HRinds)-rpp(HRinds));
    HRz(noHRinds) = 0;
    
    % H2O Transpired from each level
    TRz = UPz - (-HRz);
    TRtot = sum(UPz) - (-sum(HRz));
    
    
    %==========================================================================
    %==========================================================================
    
else    % if hydraulic redistribution is not allowed
    
    if rhc == 1     % if root hydraulic conductivity is allowed to increase
        for i = 1:nl_soil
            rootr(i) = max(0, rootfr(i)*((volliq(i)-thetadry(i))...
                /(eff_porosity(i)-thetadry(i)))); % effect of soil moisture
            % include the effect of root conductivity. Here, a linearly
            % increasing function is considered
            rootr(i) = rootr(i) + (z(i)/1000)*rootr(i);  % effect of root conductivity
        end
    else
        for i = 1:nl_soil
            rootr(i) = max(0, rootfr(i)*((volliq(i)-thetadry(i))...
                /(eff_porosity(i)-thetadry(i))));  % effect of soil moisture
        end
    end
    if (sum(rootr) > 1), rootr = rootr/sum(rootr); end
    
    % For the top soil layer
    j      = 1;
    qin    = qinfl;
    
    den    = (z(j+1)-z(j));
    num    = (smp1(j+1)-smp1(j)) - den;
    qout   = -hk(j)*num/den;
    dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den;
    dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den;
    
    rmx(j) =  qin - qout - TR*rootr(j);
    amx(j) =  0;
    bmx(j) =  dz(j)/dtime + dqodw1;
    cmx(j) =  dqodw2;
    
    q_save(1) = qin;
    
    % For the middile soil layers
    for j = 2:nl_soil - 1
        den    = (z(j) - z(j-1));
        num    = (smp1(j)-smp1(j-1)) - den;
        qin    = -hk(j-1)*num/den;
        dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den;
        dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den;
        
        den    = (z(j+1)-z(j));
        num    = (smp1(j+1)-smp1(j)) - den;
        qout   = -hk(j)*num/den;
        dqodw1 = -(-hk(j)*dsmpdw(j)   + num*dhkdw(j))/den;
        dqodw2 = -( hk(j)*dsmpdw(j+1) + num*dhkdw(j))/den;
        
        rmx(j) =  qin - qout - TR*rootr(j);
        amx(j) = -dqidw0;
        bmx(j) =  dz(j)/dtime - dqidw1 + dqodw1;
        cmx(j) =  dqodw2;
        
        q_save(j) = qin;
    end
    
    % For the bottom soil layer
    j      = nl_soil;
    den    = (z(j) - z(j-1));
    num    = (smp1(j)-smp1(j-1)) - den;
    qin    = -hk(j-1)*num/den;
    dqidw0 = -(-hk(j-1)*dsmpdw(j-1) + num*dhkdw(j-1))/den;
    dqidw1 = -( hk(j-1)*dsmpdw(j)   + num*dhkdw(j-1))/den;
    
    qout   =  hk(j);
    dqodw1 =  dhkdw(j);
    
    rmx(j) =  qin - qout - TR*rootr(j);
    amx(j) = -dqidw0;
    bmx(j) =  dz(j)/dtime - dqidw1 + dqodw1;
    cmx(j) =  0;
    
    q_save(nl_soil) = qin;
    q_save(nl_soil+1) = qout;
    
    % H2O transpired from each level
    UPz = TR*rootr;
    TRz = UPz;
    TRtot = sum(UPz);
    
    % H2O redistributed from each level
    HRz = zeros(size(UPz));
    
end

% Solve for change in soil moisture (dwat) using tridiagonal matric solver
dwat1 = TRIDIAG (nl_soil ,amx ,bmx ,cmx ,rmx);
dwat1 = dwat1(:);
% ppppp...[]

UPz = UPz(:);
HRz = HRz(:);
q_save = q_save(:);
smp1 = smp1(:);

H2O_UP = sum(UPz);
H2O_HR = sum(HRz);

% Weighted mean smp over root uptake profile [mm]
smp_weight = sum(smp1.*rootfr/sum(rootfr));% * mmH2OtoMPa;

% Weighted mean soil saturation fraction over root uptake profile
thsatfrac_weight = sum((volliq./porsl) .* UPz/sum(UPz));



dwat(:,si)=dwat1;
smp(:,si)=smp1;





