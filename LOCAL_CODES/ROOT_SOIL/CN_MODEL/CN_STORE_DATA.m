% store data


tpq=tt-NoDay_prev*(24*60/CONSTANTS.timestep);

% NUTRIENTS VARIABLES
    Cl_store(:,tpq,:)=VARIABLES.Cl;       %[gC/m3]
    Ch_store(:,tpq,:)=VARIABLES.Ch;       %[gC/m3]
    Cb_store(:,tpq,:)=VARIABLES.Cb;       %[gC/m3]
    Nl_store(:,tpq,:)=VARIABLES.Nl;        %[gN/m3] 
    Amm_store(:,tpq,:)=VARIABLES.Amm;        %[gN/m3]
    Nit_store(:,tpq,:)=VARIABLES.Nit;         %[gN/m3]
    Nitdd_store(:,tpq,:) = VARIABLES.Nitdd;                      %[g/m^3/d]deltaammdd   
    Ammdd_store(:,tpq,:) = VARIABLES.Ammdd;                      %[g/m^3/d]deltaammdd  =     
    CNl_store(:,tpq,:)=VARIABLES.CNl;        %[gC/gN]    
    dCl_dt_store(:,tpq,:) = VARIABLES.dCl_dt;  %[gC/m3/d]
    dCb_dt_store(:,tpq,:) = VARIABLES.dCb_dt;  %[gC/m3/d]
    
    kd_store(:,tpq) = PARAMS.CN.kd;
    kl_store(:,tpq) = PARAMS.CN.kl;
    kh_store(:,tpq) = PARAMS.CN.kh;
    
    
    MIN_net_store(:,tpq,:)=VARIABLES.MIN_net;     %[gN/m3/d]
    MIN_gross_store(:,tpq,:)=VARIABLES.MIN_gross;     %[gN/m3/d]
    IMM_net_store(:,tpq,:)=VARIABLES.IMM_net;      %[gN/m3/d]
    IMM_gross_store(:,tpq,:)=VARIABLES.IMM_gross;     %[gN/m3/d]
    Nreg_store(:,tpq,:)=VARIABLES.Nreg;
    DECl_store(:,tpq,:)=VARIABLES.DECl;
    PHI_store(:,tpq)=VARIABLES.PHI;
    phi_store(:,tpq)=VARIABLES.phi; 
    fSd_store(:,tpq,:)=VARIABLES.fSd;
    fTd_store(:,tpq,:)=VARIABLES.fTd;
    mberrorN_store(:,tpq)=VARIABLES.mberrorN;     %[gN/m3/d]
    mberrorC_store(:,tpq)=VARIABLES.mberrorC;     %[gC/m3/d]    
    ADD_store(:,tpq,:) = VARIABLES.ADD;            %[gr/m3/d]
    
    %new
    litterthickness_store(1,tpq,:)=VARIABLES.SOIL.litterthickness;
    
    % NITRATE
    if SWITCHES.CN.NupRootBiomass == 1
        ADV_nit_m2_store(:,tpq) = VARIABLES.ADV_nit_m2;                      %[g/m^2/d]
        ADV_nit_m3_store(:,tpq) = VARIABLES.ADV_nit_m3;                      %[g/m^3/d]
        ADVr_nit_m2_store(:,tpq) = VARIABLES.ADVr_nit_m2;                    %[g/m^2/d]
        ADVr_nit_m3_store(:,tpq) = VARIABLES.ADVr_nit_m3;                    %[g/m^3/d]
        ADVrdd_nit_m3_store(:,tpq) = VARIABLES.ADVrdd_nit_m3;                %[g/m^3/d]
        DIFF_nit_m3_store(:,tpq) = VARIABLES.DIFF_nit_m3;                    %[g/m^3/d]
        DIFF_nit_m2_store(:,tpq) = VARIABLES.DIFF_nit_m2;                    %[g/m^2/d]
    end
    LCH_nit_m3_store(:,tpq,:) = VARIABLES.LCH_nit_m3;       %[g/m^3/d]
    LCH_nit_m2_store(:,tpq,:) = VARIABLES.LCH_nit_m2;       %[g/m^2/d]
    if SWITCHES.CN.NupRootBiomass == 1
        DIFFr_nit_m2_store(:,tpq) = VARIABLES.DIFFr_nit_m2;                      %[g/m^2/d]
        DIFFr_nit_m3_store(:,tpq) = VARIABLES.DIFFr_nit_m3;                      %[g/m^3/d]
        DIFFrdd_nit_m3_store(:,tpq) = VARIABLES.DIFFrdd_nit_m3;                  %[g/m^3/d]
    end
    UP_nit_m3_store(:,tpq) = VARIABLES.UP_nit_m3;                    %[g/m^3/d]
    UP_nit_m2_store(:,tpq,:) = VARIABLES.UP_nit_m2;                    %[g/m^2/d]
    UP_nit_all_m2_store(:,tpq,:) = VARIABLES.UP_nit_all_m2;            %[g/m^2/d]
    if SWITCHES.CN.NupRootBiomass == 1
        errormbnit_store(:,tpq) = VARIABLES.errormbnit;
        deltanit_store(:,tpq) = VARIABLES.deltanit;          %[g/m^3/d]
        deltanitdd_store(:,tpq) = VARIABLES.deltanitdd;      %[g/m^3/d]
    end
    
    % Dongkook Woo - Edit
    STORAGE.UP_nit_m2_store=UP_nit_all_m2_store;
    % Dongkook Woo - Edit End

    % AMMONIUM
    if SWITCHES.CN.NupRootBiomass == 1
        ADV_amm_m2_store(:,tpq) = VARIABLES.ADV_amm_m2;                      %[g/m^2/d]
        ADV_amm_m3_store(:,tpq) = VARIABLES.ADV_amm_m3;                      %[g/m^3/d]
        ADVr_amm_m2_store(:,tpq) = VARIABLES.ADVr_amm_m2;                    %[g/m^2/d]
        ADVr_amm_m3_store(:,tpq) = VARIABLES.ADVr_amm_m3;                    %[g/m^3/d]
        ADVrdd_amm_m3_store(:,tpq) = VARIABLES.ADVrdd_amm_m3;                %[g/m^3/d]
        DIFF_amm_m3_store(:,tpq) = VARIABLES.DIFF_amm_m3;                    %[g/m^3/d]
        DIFF_amm_m2_store(:,tpq) = VARIABLES.DIFF_amm_m2;                    %[g/m^2/d]
    end
    LCH_amm_m3_store(:,tpq,:) = VARIABLES.LCH_amm_m3;       %[g/m^3/d]
    LCH_amm_m2_store(:,tpq,:) = VARIABLES.LCH_amm_m2;       %[g/m^2/d]
    if SWITCHES.CN.NupRootBiomass == 1
        DIFFr_amm_m2_store(:,tpq) = VARIABLES.DIFFr_amm_m2;                      %[g/m^2/d]
        DIFFr_amm_m3_store(:,tpq) = VARIABLES.DIFFr_amm_m3;                      %[g/m^3/d]
        DIFFrdd_amm_m3_store(:,tpq) = VARIABLES.DIFFrdd_amm_m3;                  %[g/m^3/d]
    end
    UP_amm_m3_store(:,tpq) = VARIABLES.UP_amm_m3;                    %[g/m^3/d]
    UP_amm_m2_store(:,tpq,:) = VARIABLES.UP_amm_m2;                    %[g/m^2/d]
    UP_amm_all_m2_store(:,tpq,:) = VARIABLES.UP_amm_all_m2;            %[g/m^2/d]
    if SWITCHES.CN.NupRootBiomass == 1
        errormbamm_store(:,tpq) = VARIABLES.errormbamm;
        deltaamm_store(:,tpq) = VARIABLES.deltaamm;          %[g/m^3/d]
        deltaammdd_store(:,tpq) = VARIABLES.deltaammdd;      %[g/m^3/d]
    end
    
    % Dongkook Woo - Edit
    STORAGE.UP_amm_m2_store=UP_amm_all_m2_store;
    % Dongkook Woo - Edit End
    
    % CARBON MASS BALANCE    
    Cinput_store(tpq) = VARIABLES.Cinput;
    CBoutput_store(tpq) = VARIABLES.CBoutput;
    DECout_store(tpq) = VARIABLES.DECout;
    dCb_m2_store(tpq) = VARIABLES.dCb_m2;
    dCl_m2_store(tpq) = VARIABLES.dCl_m2;    
        
    % NITROGEN MASS BALANCE   
    UP_N_m2_store(tpq,:) = VARIABLES.UP_N_m2_sp;
    LCH_N_m2_store(tpq,:) = VARIABLES.LCH_N_m2;  
    Ninput_store(tpq) = VARIABLES.Ninput;
    NBoutput_store(tpq) = VARIABLES.NBoutput;
    dNb_m2_store(tpq) = VARIABLES.dNb_m2;
    dNl_m2_store(tpq) = VARIABLES.dNl_m2;
    dAmm_m2_store(tpq,:) = VARIABLES.dAmm_m2;
    dNit_m2_store(tpq,:) = VARIABLES.dNit_m2;
    if SWITCHES.CN.NupRootBiomass == 1
        dAmmdd_m2_store(tpq) = VARIABLES.dAmmdd_m2;
        dNitdd_m2_store(tpq) = VARIABLES.dNitdd_m2;
    end
        
    % ROOT ARQUITECTURE
    if SWITCHES.CN.NupRootBiomass == 1
        RMI_store(:,tpq,:) = VARIABLES.RMI;
        RLI_store(:,tpq,:) = VARIABLES.RLI;
        RAI_store(:,tpq,:) = VARIABLES.RAI;
        RVI_store(:,tpq,:) = VARIABLES.RVI;
        RAID_store(:,tpq,:) = VARIABLES.RAID;
        RVID_store(:,tpq,:) = VARIABLES.RVID;
    end
    
    % Dongkook Woo - Edit
    % Denitrification
    Denitrif_store(:,tpq,:)  = VARIABLES.Denitrif;
    Total_N20_store(:,tpq,:) = VARIABLES.Total_N20;
    N2O_nit_store(:,tpq,:)   = VARIABLES.N2O_nit;
    
    Nitrif_store(:,tpq,:) = VARIABLES.Nitrif;
    DECh_store(:,tpq,:) = VARIABLES.DECh;
    % Dongkook Woo - End
%     
if SWITCHES.CN.Bioturbation   
    Cflux_store(tpq,:) = VARIABLES.Cbiorate;
    Nflux_store(tpq,:) = VARIABLES.Nbiorate;    
    Cl_bio_change_store(:,tpq,:) = VARIABLES.Cl_bio_change;
    Cl_bio_in_store(:,tpq,:) = VARIABLES.Cl_bio_in;   % [gr/m3/dtime]
    Cl_bio_out_store(:,tpq,:) = VARIABLES.Cl_bio_out; % [gr/m3/dtime]
    bioCerror_store(tpq,:) = VARIABLES.bioCerror;
    bioNerror_store(tpq,:) = VARIABLES.bioNerror;
    ADD_bio_store(:,tpq) = VARIABLES.ADD_bio;     %[gr/m3/d]
    ADD_ex_store(:,tpq) = VARIABLES.ADD_ex;       %[gr/m3/d]
    OUT_store(:,tpq) = VARIABLES.OUT;       % Only from bioturbation[gr/m3/d]
    ADD_net_store(:,tpq,:) = VARIABLES.ADD_net;       %[gr/m3/d]
end
%     
