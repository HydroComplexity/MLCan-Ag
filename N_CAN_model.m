
si=VARIABLES.si;

%############ Canopy Nitrogen Modeling ############%
if (tqt>CROP_GROWTH.plantingDate(si)*24 && tqt<CROP_GROWTH.harvestingDate(si)*24)
    
    VARIABLES.dim_temp(si)=VARIABLES.dim_temp(si)+1;
    
    dim_temp=VARIABLES.dim_temp(si);
    
    N_uptake_can(si) = N_uptake_can(si) + (VARIABLES.UP_N_m2_sp(si))/24;
    
    Qabs_sun_NCan(:,dim_temp,si) = VARIABLES.CANOPY.PARabs_sun;
    Tl_sun_NCan(:,dim_temp,si) = VARIABLES.CANOPY.Tl_sun(:,si);
    Ci_sun_NCan(:,dim_temp,si) = VARIABLES.CANOPY.Ci_sun(:,si);
    
    Qabs_shade_NCan(:,dim_temp,si) = VARIABLES.CANOPY.PARabs_shade;
    Tl_shade_NCan(:,dim_temp,si) = VARIABLES.CANOPY.Tl_shade(:,si);
    Ci_shade_NCan(:,dim_temp,si) = VARIABLES.CANOPY.Ci_shade(:,si);
    
    LAI_sun_can_NCan(:,dim_temp,si) = LAI_sun;
    LAI_shade_can_NCan(:,dim_temp,si) = LAI_shade;
    fLAIz_can_NCan(:,dim_temp,si) = fLAIz(:,si);
    
    Cl= CROP_GROWTH.LeafMass_can(si);
    Cr= CROP_GROWTH.RootMass_can(si);
    Cs= CROP_GROWTH.StemMass_can(si);
    Cg= CROP_GROWTH.GrainMass_can(si);
    %
    
    tqtDiff=tqt-CROP_GROWTH.plantingDate(si);
    
    %%% Here is optimization function
    % %
    if optimization==1
        if(rem(tqtDiff, time_inter*24)== 12)
%             
%             Tl_sun_NCan = Tl_sun_NCan(:,all(~isnan(Tl_sun_NCan)));
%             Tl_shade_NCan = Tl_shade_NCan(:,all(~isnan(Tl_shade_NCan)));
            
            VARIABLES_temp.NCan.Qabs_sun = mean(Qabs_sun_NCan(:,:,si),2);
            VARIABLES_temp.NCan.Tl_sun(:,si) = mean(Tl_sun_NCan(:,:,si),2);
            VARIABLES_temp.NCan.Ci_sun(:,si) = mean(Ci_sun_NCan(:,:,si),2);
            
            VARIABLES_temp.NCan.Qabs_shade = mean(Qabs_shade_NCan(:,:,si),2);
            VARIABLES_temp.NCan.Tl_shade(:,si) = mean(Tl_shade_NCan(:,:,si),2);
            VARIABLES_temp.NCan.Ci_shade(:,si) = mean(Ci_shade_NCan(:,:,si),2);
            
            VARIABLES_temp.NCan.LAI_sun_can = mean(LAI_sun_can_NCan(:,:,si),2);
            VARIABLES_temp.NCan.LAI_shade_can = mean(LAI_shade_can_NCan(:,:,si),2);
            VARIABLES_temp.NCan.fLAIz_can = mean(fLAIz_can_NCan(:,:,si),2);
            
            VARIABLES_temp.NCan.min_Vmax=VARIABLES.NCan.min_Vmax;
            VARIABLES_temp.NCan.min_chl=VARIABLES.NCan.min_chl;
            VARIABLES_temp.NCan.Nrub=VARIABLES.NCan.Nrub;
            
            VARIABLES_temp.CANOPY.Tl = VARIABLES.CANOPY.Tl;
            
            %             VARIABLES_temp.NCan.Chl=VARIABLES.NCan.Chl;
            %
            %             VARIABLES_temp.NCan.Vcmax_vert=VARIABLES.NCan.Vcmax_vert;
            %
            %             VARIABLES_temp.NCan.Nall_start=VARIABLES.NCan.Nall_start;
            
            PARAMS_temp=PARAMS;
            
            VERTSTRUC_temp=VERTSTRUC;
            
            N_uptake_can_temp=N_uptake_can(si);
            
            save ('./LOCAL_CODES/CANOPY_NITROGEN/NCAN_temp_variables.mat',...
                'PARAMS_temp', 'VARIABLES_temp', 'VERTSTRUC_temp', 'Cl',...
                'Cs','Cr', 'Cg','N_uptake_can_temp', 'CO2_Elev', 'si');
            
            MinKec=PARAMS.Nalloc.MinKec(si);
            MaxKec=PARAMS.Nalloc.MaxKec(si);
            
            MinFth=PARAMS.Nalloc.MinFth(si);
            MaxFth=PARAMS.Nalloc.MaxFth(si);
            
            MinFrub=PARAMS.Nalloc.MinFrub(si);
            MaxFrub=PARAMS.Nalloc.MaxFrub(si);
            
            
            %             cd(workpath);
            evalstr = 'Run_OptFunction';
            global Optimization_dir
            workpath = pwd;
            Optimization_dir = ['./LOCAL_CODES', '/CANOPY_NITROGEN/'];
            addpath(Optimization_dir);
            
            %% NOW RUN the model
            eval(evalstr)
            
            Qabs_sun_NCan(:,:,si) = [];
            Tl_sun_NCan(:,:,si) =  [];
            Ci_sun_NCan(:,:,si) =  [];
            
            Qabs_shade_NCan(:,:,si) = [];
            Tl_shade_NCan(:,:,si) = [];
            Ci_shade_NCan(:,:,si) =  [];
            
            LAI_sun_can_NCan(:,:,si) = [];
            LAI_shade_can_NCan(:,:,si) = [];
            fLAIz_can_NCan(:,:,si)= [];
            
            VARIABLES.dim_temp(si)=0;
            
        end
        %%%% Here end the optimization function
        
        load ('NCAN_opt_variables.mat', 'K_ec', 'K_rub_leaf', 'K_lhc_leaf');
        
    else
        %
        K_ec=0.8;
        K_rub_leaf=0.13;
        K_lhc_leaf=0.25;
        
    end
    %     disp(CROP_GROWTH.LAI_ref(tqt))
    
    VARIABLES.NCan.K_ec(si)=K_ec;
    VARIABLES.NCan.K_rub_leaf(si)=K_rub_leaf;
    VARIABLES.NCan.K_lhc_leaf(si)=K_lhc_leaf;
    
    
    LAI_ncan = VERTSTRUC.LAIz;
    
    NstrF=PARAMS.Nalloc.Nstr(si);
    NlhcF=PARAMS.Nalloc.Nlhc(si);
    Rnchl=PARAMS.Nalloc.Rnchl(si);
    
    %     LAI_ncan = LAI_ref_vert;
    
    Tl_sel=VARIABLES.CANOPY.Tl_mean(si);
    Nrub_prev=sum(VARIABLES.NCan.Nrub(:,si));
    nl_can = PARAMS.CanStruc.nl_can;
    
    [Nleaf_tot, Nr, Ns, Ng, Nleaf_la, Nrub_la, Vcmax_vert_opt,...
        Chl_opt, fn] = VCHLcal(K_ec, K_rub_leaf, K_lhc_leaf,...
        N_uptake_can(si), Nrub_prev,...
        LAI_ncan, Cl, Cr, Cs, Cg, Tl_sel, CO2_Elev, NstrF, NlhcF, Rnchl, nl_can);
    
    VARIABLES.NCan.fn(si)=fn;
    
    Q10=PARAMS.Photosyn.Q10_C4(si);
    
    Q10s = Q10.^((Tl_sel-25)./10);
    
    Vcmax_vert_opt_x = (Vcmax_vert_opt .* Q10s)./( (1+exp(0.3*(13-Tl_sel))).*(1+exp(0.3*(Tl_sel-36))) );
    
    %     Vcmax_vert_opt_x(isinf(Vcmax_vert_opt_x)|isnan(Vcmax_vert_opt_x)) = VARIABLES.NCan.min_Vmax;
    %     Vcmax_vert_opt_x(Vcmax_vert_opt_x<VARIABLES.NCan.min_Vmax)=VARIABLES.NCan.min_Vmax;
    %     Vcmax_vert_opt_x(Vcmax_vert_opt_x>VARIABLES.NCan.max_Vmax)=VARIABLES.NCan.max_Vmax;
    
    Chl_opt(isinf(Chl_opt)|isnan(Chl_opt)) = VARIABLES.NCan.min_chl;
    Chl_opt(Chl_opt<VARIABLES.NCan.min_chl)=VARIABLES.NCan.min_chl;
    Chl_opt(Chl_opt>VARIABLES.NCan.max_chl)=VARIABLES.NCan.max_chl;
    
    
    %% Save the variables
    VARIABLES.NCan.N_can(si) = N_uptake_can(si);
    VARIABLES.NCan.N_uptake_can(si) = N_uptake_can(si);
    VARIABLES.NCan.Nl(si)=Nleaf_tot;
    VARIABLES.NCan.Nr(si)=Nr;
    VARIABLES.NCan.Ns(si)=Ns;
    VARIABLES.NCan.Ng(si)=Ng;
    VARIABLES.NCan.Nleaf(:,si) = Nleaf_la;
    VARIABLES.NCan.Nrub(:,si) = Nrub_la;
    VARIABLES.NCan.Vcmax_vert(:,si)=Vcmax_vert_opt;
    
    VARIABLES.NCan.Vcmax_vert_x(:,si)=Vcmax_vert_opt_x;
    
    VARIABLES.NCan.Chl(si)=Chl_opt;
    
    VARIABLES.NCan.Nall_start(si)=1;
    
    
    %% End of Ndemand function
    
    %     K_ec_store= [K_ec_store, VARIABLES.NCan.K_ec];
    %     K_rub_leaf_store= [K_rub_leaf_store, VARIABLES.NCan.K_rub_leaf];
    %     K_lhc_leaf_store= [K_lhc_leaf_store, VARIABLES.NCan.K_lhc_leaf];
    
    
    %     fn_store=[fn_store, fn];
    
    
    %% start the Nitrogen demand
    %     rncl = 0.030;     % 0.050 (Smith et al., 2013; Johnson et al., 2007)
    %     rncr = 0.024;     % 0.040(Smith et al., 2013; Johnson et al., 2007)
    %     rncs = 0.019;     % 0.019 (Smith et al., 2013; Johnson et al., 2007)
    %     rncg = 0.019;     % 0.019 (Smith et al., 2013; Johnson et al., 2007)
    
    rncl = 0.050;       %(Smith et al., 2013; Johnson et al., 2007)
    rncr = 0.040;       %(Smith et al., 2013; Johnson et al., 2007)
    rncs = 0.019;       %(Smith et al., 2013; Johnson et al., 2007)
    rncg = 0.019;       %(Smith et al., 2013; Johnson et al., 2007)
    
    N_Demand = ((CROP_GROWTH.NPPL(si) * rncl) + (CROP_GROWTH.NPPS(si) * rncs) +...
        (CROP_GROWTH.NPPR(si) * rncr) + (CROP_GROWTH.NPPG(si) * rncg)); % gm/m2/hour
    
    N_Demand(isinf(N_Demand)|isnan(N_Demand)) = 0;
    N_Demand(N_Demand<0) = 0;
    
    CROP_GROWTH.N_Demand(si)=N_Demand;
    
else
    
    N_uptake_can(si) =0;
    
    VARIABLES.NCan.N_can(si) = 0.0;
    VARIABLES.NCan.N_uptake_can(si) =0;
    VARIABLES.NCan.Nl(si)=0;
    VARIABLES.NCan.Nr(si)=0;
    VARIABLES.NCan.Ns(si)=0;
    VARIABLES.NCan.Ng(si)=0;
    VARIABLES.NCan.Nleaf(:,si) = zeros (nl_can, 1);
    VARIABLES.NCan.Nrub(:,si) = zeros (nl_can, 1);
    VARIABLES.NCan.Vcmax_vert(:,si) = zeros (nl_can, 1);
    VARIABLES.NCan.Vcmax_vert_x(:,si) = zeros (nl_can, 1);
    
    VARIABLES.NCan.Chl(si) = 0;
    
    VARIABLES.NCan.Nall_start(si)=0;
    
    CROP_GROWTH.N_Demand(si)=0;
    
    VARIABLES.NCan.K_ec(si)=0;
    VARIABLES.NCan.K_rub_leaf(si)=0;
    VARIABLES.NCan.K_lhc_leaf(si)=0;
    VARIABLES.NCan.fn(si)=0;
    
end

if (tqt>CROP_GROWTH.plantingDate(si)*24 && tqt<CROP_GROWTH.doy3(si)*24)
    if(sum(LAI_ncan)<0.5)
        VARIABLES.NCan.Nall_start(si)=0;
    end
end

