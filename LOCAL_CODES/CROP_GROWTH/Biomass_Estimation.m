function [CROP_GROWTH ] = Biomass_Estimation( CROP_GROWTH, CONSTANTS, PARAMS, Ph_can, Tl, doy1, doy2, doy3, doy4, doy5, si, t)

%
%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%%                              FUNCTION CODE                            %%
%%         Estimation of biomass and LAI using allometric approach       %%


%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Rohit Nandan                                           %
%   Date         : November 18, 2019                                      %
%-------------------------------------------------------------------------%

%BIOMASS_ESTIMATION Summary of this function goes here
%   Detailed explanation goes here

WSTRES = 0;

plantingDate = CROP_GROWTH.plantingDate(si);
harvestingDate = CROP_GROWTH.harvestingDate(si);

LEAFMASS_can = CROP_GROWTH.LeafMass_can(si);
STEMMASS_can = CROP_GROWTH.StemMass_can(si);
ROOTMASS_can = CROP_GROWTH.RootMass_can(si);
GRAINMASS_can = CROP_GROWTH.GrainMass_can(si);

LAIsim_can = CROP_GROWTH.LAIsim_can(si);

CF_LF=PARAMS.CGM.CF_LF(si);
CF_ST=PARAMS.CGM.CF_ST(si);
CF_RT=PARAMS.CGM.CF_RT(si);
CF_GN=PARAMS.CGM.CF_GN(si);

HI=PARAMS.CGM.HI(si);

uca = 10.^-6;
ucb = 10.^-3;

td = 24 * 60 / CONSTANTS.timestep;

tq = t + (CROP_GROWTH.DOY_start-1) * td;

incr_factor = 1;

PSN_can = Ph_can;                          % +Photosynthesis;

%      disp(PSN_can);

%   IPA = SWITCHES.IPA;
IPA = 1;

DILE_FC = 0;
DILE_FW = 0;
LF_OVRC = 0;
ST_OVRC = 0;
RT_OVRC = 0;
LFPT = 0;
STPT = 0;
RTPT = 0;
GRAINPT = 0;

CBHydraFx_can = PSN_can * incr_factor * 30 * (1 - WSTRES) * IPA * uca;


if tq >= (doy1-1) * td && tq < (doy2-1) * td				%// PGS 3
    LFPT = PARAMS.CGM.LFPT_3(si);
    STPT = PARAMS.CGM.STPT_3(si);
    RTPT = PARAMS.CGM.RTPT_3(si);
    
elseif tq >= (doy2-1) * td && tq < (doy3-1) * td			%// PGS 4
    LFPT = PARAMS.CGM.LFPT_4(si);
    STPT = PARAMS.CGM.STPT_4(si);
    RTPT = PARAMS.CGM.RTPT_4(si);
    
elseif tq >= (doy3-1) * td && tq < (doy4-1) * td			%// PGS 5
    DILE_FC = PARAMS.CGM.DILE_FC_5(si);
    DILE_FW = PARAMS.CGM.DILE_FW_5;
    LF_OVRC = PARAMS.CGM.LF_OVRC_5(si);
    ST_OVRC = PARAMS.CGM.ST_OVRC_5(si);
    RT_OVRC = PARAMS.CGM.RT_OVRC_5(si);
    
    LFPT = PARAMS.CGM.LFPT_5(si);
    RTPT = PARAMS.CGM.RTPT_5(si);
    STPT = PARAMS.CGM.STPT_5(si);
    GRAINPT = PARAMS.CGM.GRAINPT_5(si);
    
elseif tq >= (doy4-1) * td && tq < (doy5-1) * td			%// PGS 6
    DILE_FC = PARAMS.CGM.DILE_FC_6(si);
    DILE_FW = PARAMS.CGM.DILE_FW_6;
    LF_OVRC = PARAMS.CGM.LF_OVRC_6(si);
    ST_OVRC = PARAMS.CGM.ST_OVRC_6(si);
    RT_OVRC = PARAMS.CGM.RT_OVRC_6(si);
    
    LFPT = PARAMS.CGM.LFPT_6(si);
    RTPT = PARAMS.CGM.RTPT_6(si);
    STPT = PARAMS.CGM.STPT_6(si);
    GRAINPT = PARAMS.CGM.GRAINPT_6(si);
    
end

% 	//Maintenance Respiration///////////////////////
FNF = 1;

TF_can = PARAMS.CGM.Q10MR(si).^ ((Tl - 25.16) / 10);

% 	RESP_can = CROP_GROWTH.LFMR25 * FNF * TF_can * (CROP_GROWTH.LAI_actual)*(1 - WSTRES);
RESP_can = PARAMS.CGM.LFMR25(si) * FNF * TF_can * (LAIsim_can)*(1 - WSTRES);
RSLEAF_can = RESP_can * 30 * uca;
RSSTEM_can = PARAMS.CGM.STMR25(si)*(STEMMASS_can*ucb)*TF_can * 30 * uca;
RSROOT_can = PARAMS.CGM.RTMR25(si)*(ROOTMASS_can*ucb)*TF_can * 30 * uca;
RSGRAIN_can = PARAMS.CGM.GRAINMR25(si)*(GRAINMASS_can*ucb)*TF_can * 30 * uca;

GRLEAF_can = PARAMS.CGM.FRA_GR(si) * ((LFPT * CBHydraFx_can) - RSLEAF_can);
GRSTEM_can = PARAMS.CGM.FRA_GR(si) * ((STPT * CBHydraFx_can) - RSSTEM_can);
GRROOT_can = PARAMS.CGM.FRA_GR(si) * ((RTPT * CBHydraFx_can) - RSROOT_can);
GRGRAIN_can = PARAMS.CGM.FRA_GR(si) * ((GRAINPT * CBHydraFx_can) - RSGRAIN_can);



% 	//Turnover and Death//////////////////

LFTOVR_can = LF_OVRC * LEAFMASS_can * uca;
STTOVR_can = ST_OVRC * STEMMASS_can * uca;
RTTOVR_can = RT_OVRC * ROOTMASS_can * uca;


SC_can = exp(-0.3*((Tl + 273) - PARAMS.CGM.LEFREEZ(si)))*(LEAFMASS_can / 120);
% 	//SC_can = 0;
SD_can = exp((WSTRES - 1) * (1 - WSTRES));
% 	//SD_can = 0;
DIELF_can = (LEAFMASS_can) * uca*(DILE_FW * SD_can + DILE_FC * SC_can);


% 	// Carbohydrate flux allocation/////////////////

NPPL_can = (LFPT * CBHydraFx_can)- (GRLEAF_can)-(RSLEAF_can);
NPPS_can = (STPT * CBHydraFx_can)- (GRSTEM_can)-(RSSTEM_can);
NPPR_can = (RTPT * CBHydraFx_can)- (GRROOT_can)-(RSROOT_can);
NPPG_can = (GRAINPT * CBHydraFx_can)- (GRGRAIN_can)-(RSGRAIN_can);


%     disp(NPPL_can);

%     if NPPL_can<0
%         NPPL_can = 0;
%     end
%
% 	if NPPS_can<0
%         NPPS_can = 0;
%     end
%
% 	if NPPR_can<0
%         NPPR_can = 0;
%     end
%
% 	if NPPG_can<0
%         NPPG_can = 0;
%     end


% 	// Total Biomass accumulation///////////////////////

% 	//LEAFMASS_can = LEAFMASS_can + (NPPL_can - LFTOVR_can) * CONSTANTS->dtime;						// Assumption that leaf will not be affected by water and heat stress

dtime=CONSTANTS.dtime;

LEAFMASS_can = LEAFMASS_can + (NPPL_can - LFTOVR_can - DIELF_can) * dtime;
STEMMASS_can = STEMMASS_can + (NPPS_can - STTOVR_can) * dtime;
ROOTMASS_can = ROOTMASS_can + (NPPR_can - RTTOVR_can) * dtime;
GRAINMASS_can = GRAINMASS_can + (NPPG_can)* dtime;

if tq > plantingDate*td && tq < harvestingDate*td
    CROP_GROWTH.LeafMass_can(si) = LEAFMASS_can;
    CROP_GROWTH.StemMass_can(si) = STEMMASS_can;
    CROP_GROWTH.RootMass_can(si) = ROOTMASS_can;
    CROP_GROWTH.GrainMass_can(si) = GRAINMASS_can;
    
else
    
    CROP_GROWTH.LeafMass_can(si) = 0.0;
    CROP_GROWTH.StemMass_can(si) = 0.0;
    CROP_GROWTH.RootMass_can(si) = 0.0;
    CROP_GROWTH.GrainMass_can(si) = 0.0;
    CROP_GROWTH.LAIsim_can(si) = 0.0;
    
end

if isnan(CROP_GROWTH.LeafMass_can(si)) | CROP_GROWTH.LeafMass_can(si) < 0.0
    CROP_GROWTH.LeafMass_can(si) = 0.0;
end
if isnan(CROP_GROWTH.StemMass_can(si)) | CROP_GROWTH.StemMass_can(si) < 0.0
    CROP_GROWTH.StemMass_can(si) = 0.0;
end
if isnan(CROP_GROWTH.GrainMass_can(si)) | CROP_GROWTH.GrainMass_can(si) < 0.0
    CROP_GROWTH.GrainMass_can(si) = 0.0;
end
if isnan(CROP_GROWTH.RootMass_can(si)) | CROP_GROWTH.RootMass_can(si) < 0.0
    CROP_GROWTH.RootMass_can(si) = 0.0;
end


AboveMass_can(si) = CROP_GROWTH.LeafMass_can(si) + CROP_GROWTH.StemMass_can(si) + CROP_GROWTH.GrainMass_can(si);
TotalMass_can(si) = AboveMass_can(si) + CROP_GROWTH.RootMass_can(si);


CROP_GROWTH.LeafBiomass(si) = CROP_GROWTH.LeafMass_can(si)/CF_LF;
CROP_GROWTH.StemBiomass(si) = CROP_GROWTH.StemMass_can(si)/CF_ST;
CROP_GROWTH.RootBiomass(si) = CROP_GROWTH.RootMass_can(si)/CF_RT;
CROP_GROWTH.GrainBiomass(si) = CROP_GROWTH.GrainMass_can(si)/CF_GN;

CROP_GROWTH.AboveBiomass_can(si) = CROP_GROWTH.LeafBiomass(si) + CROP_GROWTH.StemBiomass(si) + CROP_GROWTH.GrainBiomass(si);
CROP_GROWTH.TotalBiomass_can(si) = CROP_GROWTH.AboveBiomass_can(si) + CROP_GROWTH.RootMass_can(si);

CROP_GROWTH.GrainYield(si)=HI * CROP_GROWTH.GrainBiomass(si);

LAIsim_can = CROP_GROWTH.LeafBiomass(si) * PARAMS.CGM.BIO2LAI(si);

LAIsim_can(LAIsim_can<=0.1)=0.1;

LAIsim_can(LAIsim_can>6.5)=6.5;


CROP_GROWTH.LAIsim_can(si) = LAIsim_can;

CROP_GROWTH.DeathRate(si) = (STTOVR_can + LFTOVR_can + DIELF_can)  * dtime;

% disp(NPPL_can);

CROP_GROWTH.NPPL(si) = NPPL_can * dtime;
CROP_GROWTH.NPPS(si) = NPPS_can * dtime;
CROP_GROWTH.NPPR(si) = NPPR_can * dtime;
CROP_GROWTH.NPPG(si) = NPPG_can * dtime;

CROP_GROWTH.RootDeathRate(si) = RTTOVR_can  * dtime;


end


