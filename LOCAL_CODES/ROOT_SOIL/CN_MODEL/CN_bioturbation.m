
function [Cl, VARIABLES] = CN_bioturbation (PARAMS, VARIABLES, CONSTANTS,...
    FORCING, VERTSTRUC, SWITCHES, fTd, si)

Cl_2=VARIABLES.Cl(:,si);

% ALLOCATE MATRICES TO USE
nspecies = PARAMS.CanStruc.nspecies;
% Dongkook Woo - Edit
if SWITCHES.CN.NupRootBiomass == 1
    CNveg = nan;
elseif SWITCHES.CN.NupRootBiomass == 0
    CNabove = nan;
    CNbelow = nan;
end
% Dongkook Woo - Edit End

%   INPUTS:
% DE REFERENCE BLOCKS
%  VARIABLES structure
timestep = VARIABLES.timestep;          % timestep = Current time step
DECl = VARIABLES.DECl;                  % Decomposition in litter layer
BD = VARIABLES.BD;                      % Biomass death    
litter_dz = VARIABLES.SOIL.litterthickness(si);%   layer thickness [m]
CNl = VARIABLES.CNl(:,si);
%  FORCING structure
% Dongkook Woo - Edit
if SWITCHES.CN.NupRootBiomass == 1
%     TBMla = FORCING.TBMla(:,timestep);      % Input of litter in current time step for each species [gr/m2]
    TBMla = VARIABLES.CANOPY.DeathRate(si);
elseif SWITCHES.CN.NupRootBiomass == 0
%     TBMla = FORCING.BMla(:,timestep);      % Input of litter in current time step for each species [gr/m2]
    TBMla = VARIABLES.CANOPY.DeathRate(si);
end
% Dongkook Woo - Edit End
ADDa = sum(TBMla);                      % Total input of litter in current time step [gr / m^2]  
%  PARAMS structure     
% for ii=1:1:nspecies
% Dongkook Woo - Edit
if SWITCHES.CN.NupRootBiomass == 1
    CNveg = PARAMS.CN.CNveg{si}(timestep);
elseif SWITCHES.CN.NupRootBiomass == 0
    CNabove = FORCING.CNabove(si,timestep)';
end
% Dongkook Woo - Edit End
% end

% Dongkook Woo - Edit
if SWITCHES.CN.NupRootBiomass == 1
    CNveg_mean = sum(CNveg'.*TBMla/sum(TBMla));
elseif SWITCHES.CN.NupRootBiomass == 0
    CNveg_mean = sum(CNabove'.*TBMla/sum(TBMla));
    if isnan(CNveg_mean)
        CNveg_mean = mean(CNabove);
    end
end

% Dongkook Woo - Edit End
kbio = PARAMS.kbio;
nl_soil=PARAMS.nl_soil;                 % Number of layers in the soil
Clitter = PARAMS.CN.Clitter(si);
% Dongkook - Comment
%layerbio = PARAMS.CN.layerbio;
% Dongkook - Comment End
Dtop= PARAMS.CN.D;
fHorO = PARAMS.CN.fHorO;
% Dongkook Woo - Edit
%fTdl = PARAMS.CN.fTdl(timestep);        % Temperature factor in litter.
fTdl = VARIABLES.fTdl(si);
% Dongkook Woo - Edit End
%  VERTSTRUC structure
dz = VERTSTRUC.dzs;
zh = VERTSTRUC.zhs;                     % nodes depth in the grid [m]   
zn = VERTSTRUC.zns;                     % nodes depth in the grid [m]   
%  CONSTANTS structure
dtime = CONSTANTS.dtime;                % dtime = Time step [1800 s]

% BIOTURBATION AND LITTER FRAGMENTATION FLUX FROM LITTER TO O HORIZON 
%**************************************************************************
%PotBiot = PotBiot/(365*3600*24)*dtime;   %[grC/m2] in the time step 
Biotflux = Clitter*litter_dz*kbio*fTdl; %[grC/m2/dtime] in the time step 
Cflux = min(Biotflux,Cl_2(1)*(litter_dz-0.005));     %[grC/m2/dtime] in the time step 

% change in total carbon in litter per m2
decom = (DECl-BD)*litter_dz*dtime/86400;  % [grC/m2/dtime]
%decom = (Clitter-Cl(1))*litter_dz;    % just to check both should be the same

dClit_m2 = ADDa - Cflux - decom;    %[grC/m2]

% Conmpute new CN ratio in litter layer
CNl(1) = (Cl_2(1)*litter_dz*CNl(1) + dClit_m2(1)*CNveg_mean)/(Cl_2(1)*litter_dz + dClit_m2(1));

% change in litter thickness
litterthickness = (Cl_2(1)*litter_dz + dClit_m2(1))/Clitter;  % [m] Satisfy constant litter concentration
Cl_2(1) = Clitter;                                          % Set litter concentration
deltathick = litterthickness - VARIABLES.SOIL.litterthickness(si); % [m]
VARIABLES.SOIL.litterthickness(si) = litterthickness;    % [m]


%               BIOTURBATION DUE TO DIFFUSION 
%**************************************************************************

BC = zeros(nl_soil,1);
BC(1) = Cflux*(1-fHorO);       % Bioturbation rate 1st layer [grC/m2/dtime]
BC(2) = Cflux*fHorO;       % Bioturbation rate 2nd layer [grC/m2/dtime]
BC(nl_soil) = 0;

% Bioturbation at the top surface
D = Dtop*exp(-0.1*(zh*100));      % [cm2/year]  
% Compute Dbio at each interface layer exponential decrease function from Cousins etal Chemosphere 1999 
% Dongkook Woo - Edit
%D = D/100/100/365/48;          % from [cm2/year] to [m2/dtime]
D = D/100/100/365/(24*60*60/dtime);          % from [cm2/year] to [m2/dtime]
% Dongkook Woo - Edit End

%D(3:12) =0;
% Compute distance between nodes
deltaz = zn(2:nl_soil) - zn(1:nl_soil-1);

%Compute the matrices
[A, B]  = CN_biomatrices  (dz, deltaz, D, 1, BC);

% 
% Compute new states of C
Clnew = A^(-1)*(Cl_2(2:nl_soil+1)+B);
% Clnew = Clnew(:);
% pppppp;
% Compute fluxes in and fluxes out in each layer
[Cin_m2, Cout_m2, difbio_m2, Bioflux] = CN_biofluxes (Clnew, dz, deltaz, D, BC);

% Expand the matrixes to influde litter change to dauly
Bioflux = [ [0 0]; Bioflux]*86400/dtime;  %[gr/m2/d]
Cin_m2 = [ 0 ; Cin_m2]*86400/dtime;   %[gr/m2/d]
Cout_m2 = [ 0 ; Cout_m2]*86400/dtime; %[gr/m2/d]
difbio_m2 = [ 0 ; difbio_m2]*86400/dtime;  %[gr/m2/d]
 
% Change units to [gr/m3/d]
difbio_m3 = difbio_m2./[litterthickness ; dz];  %[gr/m3/d]
Cin_m3 = Cin_m2./[litterthickness ; dz];       %[gr/m3/d]
Cout_m3 = Cout_m2./[litterthickness ; dz];     %[gr/m3/d]

VARIABLES.Cl_bio_change(:,si) = difbio_m3;                % [gr/m3/d]
VARIABLES.Cl_bio_in(:,si) = Cin_m3;                        % [gr/m3/d]                    
VARIABLES.Cl_bio_out(:,si) = Cout_m3;                      % [gr/m3/d]

% Create weighted fraction of the positive values coming into each layer 
indmatriz = zeros(length(Cin_m3),2);
indmatriz(Bioflux>0) = Bioflux(Bioflux>0);
wgtCinput = (indmatriz./repmat(sum(indmatriz,2),1,2));
indnan = isnan(sum(wgtCinput,2));
% CN 
matrixinCN = [[CNveg_mean ; CNl(1:nl_soil)] [CNl(2:nl_soil+1) ; CNl(nl_soil+1)]];
CN_bio_in = (matrixinCN(:,1).*matrixinCN(:,2))./((matrixinCN(:,2).*wgtCinput(:,1))+((matrixinCN(:,1).*wgtCinput(:,2))));
%CN_bio_in = sum(matrixinCN.*wgtCinput,2);
%CN_bio_in = (wgtCinput(:,1).*Cin_m2(:)+wgtCinput(:,2).*Cin_m2(:))./...
%    (wgtCinput(:,1).*Cin_m2(:)./matrixinCN(:,1)+wgtCinput(:,2).*Cin_m2(:)./matrixinCN(:,2));
%N_in = (wgtCinput(:,1).*Cin_m2(:)./matrixinCN(:,1)+wgtCinput(:,2).*Cin_m2(:)./matrixinCN(:,2));
%N_out = 
CN_bio_in(indnan) = CNl(indnan);
% CN
CN_bio_out = CNl;
VARIABLES.CN_bio_in(:,si) = CN_bio_in;
VARIABLES.CN_bio_out(:,si) = CN_bio_out;

%**************************************************************************
% CHECK MASS BALANCE BY SOLUTION OF BIOTURBATION

% CARBON
cinput = Cflux.*86400/dtime;     % [grN/m2/d]
Cchange = (Cin_m3-Cout_m3).*[litterthickness ; dz]; %[grN/m2/d]
VARIABLES.bioCerror(:,si) = (cinput-sum(Cchange));  %[grC/m2/d] 

% NITROGEN
ninput = Cflux/(CNl(1)).*86400/dtime;    % [grN/m2/d]
Nin = (Cin_m3.*[litterthickness ; dz])./CN_bio_in;
Nout = (Cout_m3.*[litterthickness ; dz])./CN_bio_out;
Nchange = (Cin_m3./CN_bio_in - Cout_m3./CN_bio_out).*[litterthickness ; dz]; %[grN/m2/dtime]
VARIABLES.bioNerror(:,si) = (ninput-sum(Nchange))*86400/dtime;  %[grN/m2/d] 

% SAVE Bioturbation and fragmentation rateS from litter to Horizon.
VARIABLES.Cbiorate(:,si) = cinput;    
VARIABLES.Nbiorate(:,si) = ninput;    

Cl(:,si)=Cl_2;


