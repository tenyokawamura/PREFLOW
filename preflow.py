from preflow_h import *
from scipy import integrate

def preflow(engs, params, fluxes):
    # ------------------------------------- #
    # ---------- Initial setting ---------- #
    # ------------------------------------- #
    inpar=SetParameter()
    inpar.set_inpar(pars=params, es=engs)
    
    #######################
    ### Energy spectrum ###
    #######################
    if inpar.quant==0:
        specx=SpectralModelXspec()
        fluxes_d =specx.diskbb_spec   (es=inpar.es, temp=inpar.temp_bb_d, norm=inpar.norm_d)
        fluxes_s =specx.nthcomp_spec  (es=inpar.es, gamma=inpar.gamma_s, ktbb=inpar.temp_bb_s, kte=inpar.temp_e_s, norm=inpar.norm_s)
        fluxes_h =specx.nthcomp_spec  (es=inpar.es, gamma=inpar.gamma_h, ktbb=inpar.temp_bb_h, kte=inpar.temp_e_h, norm=inpar.norm_h)
        fluxes_sr=specx.relxillcp_spec(\
            es=inpar.es,         incl=inpar.incl,      a=inpar.a,\
            r_in=inpar.r_in_sr,  r_out=inpar.r_out_sr, index=inpar.index_sr,\
            gamma=inpar.gamma_s, logxi=inpar.logxi_sr, logcn=inpar.logcn_sr,\
            cafe=inpar.cafe,     kte=inpar.temp_e_s,   cfref=inpar.cf_sr)
        fluxes_hr=specx.relxillcp_spec(\
            es=inpar.es,         incl=inpar.incl,      a=inpar.a,\
            r_in=inpar.r_in_hr,  r_out=inpar.r_out_hr, index=inpar.index_hr,\
            gamma=inpar.gamma_h, logxi=inpar.logxi_hr, logcn=inpar.logcn_hr,\
            cafe=inpar.cafe,     kte=inpar.temp_e_h,   cfref=inpar.cf_hr)

        n_e=len(inpar.es)
        for i_e in range(n_e-1):
            fluxes[i_e]=fluxes_d[i_e]+fluxes_s[i_e]+fluxes_h[i_e]+fluxes_sr[i_e]+fluxes_hr[i_e]

    #########################
    ### Timing properties ###
    #########################
    elif inpar.quant in [1, 2, 3, 4, 5, 6]:
        # ------------------------------------ #
        # ---------- Set basic unit ---------- #
        # ------------------------------------ #
        bunit=BasicUnit()
        bunit.set_unit(mass=inpar.mass)
        # ----- Print basic information ----- #
        if inpar.display==1:
            print('--------------------------------------------------')
            print('M   : {0:.1f} solar mass'.format(bunit.m))
            print('Rg  : {0:.1f} km'.format(bunit.rg))

        # ---------------------------------------- #
        # ---------- Set time/frequency ---------- #
        # ---------------------------------------- #
        tifr=TimeFrequency()
        f_data_min=inpar.fs_data[0]
        f_data_max=inpar.fs_data[-1]
        tifr.set_par(f_data_min=f_data_min, f_data_max=f_data_max)
        tifr.t_set()
        tifr.f_set()

        # ----- Print timing information ----- #
        if inpar.display==1:
            print('--------------------------------------------------')
            print('Timing interval: {0:.0e} sec'.format(tifr.dt))
            print('Number of data: {0:.0f}'.format(tifr.n_data))

        # ----------------------------------------------- #
        # ---------- Accretion flow separation ---------- #
        # ----------------------------------------------- #
        rings=Flow2Ring()
        rings.n_ring=int(inpar.n_ring)
        rings.r_in =inpar.r_in
        rings.r_out=inpar.r_out
        rings.par_set()
        rings.coeff_calc()
        rings.lowboundr_calc()
        rings.centerr_calc()
        rings.width_calc()
        rings.interval_calc()
        rings.n_ring_dec_calc()

        # ----------------------------------------------------- #
        # ---------- Assign properties to each ring  ---------- #
        # ----------------------------------------------------- #
        r_min_eps=inpar.r_in
        # --- Viscous frequency --- #
        is_flow =np.where(rings.rs_min<=inpar.r_ds)[0]
        is_hcomp=np.where(rings.rs_min<=inpar.r_sh)[0]
        is_scomp=np.where((inpar.r_sh<rings.rs_min) & (rings.rs_min<=inpar.r_ds))[0]
        is_disk =np.where(inpar.r_ds<rings.rs_min)[0]
        rs_flow =rings.rs[is_flow]
        rs_hcomp=rings.rs[is_hcomp]
        rs_scomp=rings.rs[is_scomp]
        rs_disk =rings.rs[is_disk]
        # Width is used to calculate weight. (2022/01/28)
        wids_flow =rings.wids[is_flow]
        wids_hcomp=rings.wids[is_hcomp]
        wids_scomp=rings.wids[is_scomp]
        wids_disk =rings.wids[is_disk]

        # Variability time-scale
        fs_vis_flow=f_vis_calc(r=rs_flow, lb=inpar.lb_flow, m=inpar.m_flow) #[c/Rg]
        fs_vis_disk=f_vis_calc(r=rs_disk, lb=inpar.lb_disk, m=inpar.m_disk) #[c/Rg]
        rings.fs_vis=np.append(fs_vis_flow, fs_vis_disk)
        # Propagation time-scale
        fs_acc_flow=f_vis_calc(r=rs_flow, lb=inpar.lb_acc_flow, m=inpar.m_acc_flow) #[c/Rg]
        fs_acc_disk=f_vis_calc(r=rs_disk, lb=inpar.lb_acc_disk, m=inpar.m_acc_disk) #[c/Rg]
        rings.fs_acc=np.append(fs_acc_flow, fs_acc_disk)
        rings.vs_rad=rings.rs*rings.fs_acc #[c]

        ### Unecessary anymore, should be removed in the future ###
        # --- Damping factor --- # 
        cds_flow=(inpar.cd_flow**(1./rings.n_dec))*np.ones(len(rs_flow))
        cds_disk=(inpar.cd_disk**(1./rings.n_dec))*np.ones(len(rs_disk))
        if len(rs_disk)==0:
            cds_flow[-1]=1.
        else:
            cds_flow[-1]=inpar.cd_tran
            cds_disk[-1]=1.
        rings.cds=np.append(cds_flow, cds_disk)

        # --- Variability amplitude --- #
        lfs_var_flow=gauss(x=rs_flow, norm=inpar.lf_var_flow, mu=inpar.r_ds, sigma=inpar.r_sig_flow)
        lfs_var_disk=gauss(x=rs_disk, norm=inpar.lf_var_disk, mu=inpar.r_ds, sigma=inpar.r_sig_disk)
        rings.lfs_var=np.append(lfs_var_flow, lfs_var_disk)

        # --- Final operation --- #
        rings.out2in()
        rs_flow =rs_flow[::-1]
        rs_hcomp=rs_hcomp[::-1]
        rs_scomp=rs_scomp[::-1]
        rs_disk =rs_disk[::-1]
        wids_flow =wids_flow[::-1]
        wids_hcomp=wids_hcomp[::-1]
        wids_scomp=wids_scomp[::-1]
        wids_disk =wids_disk[::-1]

        # --- Weight --- #
        # Stressed
        if inpar.stress==1:
            stress=True 
        # Stress-free
        elif inpar.stress==2:
            stress=False
        else:
            print('Error')
            sys.exit()

        ### Energy band ###
        # Spectrum
        specx=SpectralModelXspec()
        es=np.array([inpar.e_min, inpar.e_max])
        de=inpar.e_max-inpar.e_min
        flux_d =specx.diskbb_spec (es=es, temp=inpar.temp_bb_d, norm=inpar.norm_d)[0]/de
        flux_s =specx.nthcomp_spec(es=es, gamma=inpar.gamma_s, ktbb=inpar.temp_bb_s, kte=inpar.temp_e_s, norm=inpar.norm_s)[0]/de
        flux_h =specx.nthcomp_spec(es=es, gamma=inpar.gamma_h, ktbb=inpar.temp_bb_h, kte=inpar.temp_e_h, norm=inpar.norm_h)[0]/de
        fluxr_d=0.
        fluxr_s=specx.relxillcp_spec(\
            es=es,               incl=inpar.incl,      a=inpar.a,\
            r_in=inpar.r_in_sr,  r_out=inpar.r_out_sr, index=inpar.index_sr,\
            gamma=inpar.gamma_s, logxi=inpar.logxi_sr, logcn=inpar.logcn_sr,\
            cafe=inpar.cafe,     kte=inpar.temp_e_s,   cfref=inpar.cf_sr)[0]/de
        fluxr_h=specx.relxillcp_spec(\
            es=es,               incl=inpar.incl,      a=inpar.a,\
            r_in=inpar.r_in_hr,  r_out=inpar.r_out_hr, index=inpar.index_hr,\
            gamma=inpar.gamma_h, logxi=inpar.logxi_hr, logcn=inpar.logcn_hr,\
            cafe=inpar.cafe,     kte=inpar.temp_e_h,   cfref=inpar.cf_hr)[0]/de
        #print(flux_d, flux_s, flux_h, fluxr_d, fluxr_s, fluxr_h)

        # Normalize such that total corresponds to unity.
        flux_tot=flux_d+flux_s+flux_h+fluxr_d+fluxr_s+fluxr_h
        flux_d /=flux_tot
        flux_s /=flux_tot
        flux_h /=flux_tot
        fluxr_d/=flux_tot
        fluxr_s/=flux_tot
        fluxr_h/=flux_tot
        #print(flux_d, flux_s, flux_h, fluxr_d, fluxr_s, fluxr_h)

        # eta (average)
        # No upper and lower bounds
        eta_d =calc_eta_ave(e_min=inpar.e_min, e_max=inpar.e_max, eta0=inpar.eta0_d, eta1=inpar.eta1_d)
        eta_s =calc_eta_ave(e_min=inpar.e_min, e_max=inpar.e_max, eta0=inpar.eta0_sc, eta1=inpar.eta1_sc)
        eta_h =calc_eta_ave(e_min=inpar.e_min, e_max=inpar.e_max, eta0=inpar.eta0_hc, eta1=inpar.eta1_hc)
        etar_d=0.
        etar_s=calc_eta_ave(e_min=inpar.e_min, e_max=inpar.e_max, eta0=inpar.eta0_sr, eta1=inpar.eta1_sr)
        etar_h=calc_eta_ave(e_min=inpar.e_min, e_max=inpar.e_max, eta0=inpar.eta0_hr, eta1=inpar.eta1_hr)
        #print(eta_d, eta_s, eta_h, etar_d, etar_s, etar_h)

        # \eta*C(E)
        cc_d =flux_d *eta_d
        cc_s =flux_s *eta_s
        cc_h =flux_h *eta_h
        ccr_d=fluxr_d*etar_d
        ccr_s=fluxr_s*etar_s
        ccr_h=fluxr_h*etar_h

        # Energy dissipated
        es_dis=epsilon_calc(r=rs_hcomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
            2.*np.pi*rs_hcomp*wids_hcomp
        # Total energy dissipated (normalization)
        es_dis_tot=np.sum(es_dis)
        # Direct
        ws_hcomp=cc_h*es_dis/es_dis_tot
        # Reprocessed
        wrs_hcomp=ccr_h*es_dis/es_dis_tot

        es_dis=epsilon_calc(r=rs_scomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
            2.*np.pi*rs_scomp*wids_scomp
        es_dis_tot=np.sum(es_dis)
        # Direct
        ws_scomp=cc_s*es_dis/es_dis_tot
        # Reprocessed
        wrs_scomp=ccr_s*es_dis/es_dis_tot

        es_dis=epsilon_calc(r=rs_disk, stress=stress,  gamma=inpar.gamma_disk, r_min=inpar.r_min)*\
            2.*np.pi*rs_disk*wids_disk
        es_dis_tot=np.sum(es_dis)
        # Direct
        ws_disk=cc_d*es_dis/es_dis_tot
        # Reprocessed
        wrs_disk=ccr_d*es_dis/es_dis_tot

        # Direct
        ws=ws_disk
        ws=np.append(ws, ws_scomp)
        ws=np.append(ws, ws_hcomp)
        # Reprocessed
        wrs=wrs_disk
        wrs=np.append(wrs, wrs_scomp)
        wrs=np.append(wrs, wrs_hcomp)

        #print(ws_disk)
        #print(wrs_disk)
        #print(ws)
        #print(wrs)

        ### Reference band ###
        # Spectrum
        specx=SpectralModelXspec()
        es=np.array([inpar.e_minr, inpar.e_maxr])
        de=inpar.e_maxr-inpar.e_minr
        flux_d =specx.diskbb_spec (es=es, temp=inpar.temp_bb_d, norm=inpar.norm_d)[0]/de
        flux_s =specx.nthcomp_spec(es=es, gamma=inpar.gamma_s, ktbb=inpar.temp_bb_s, kte=inpar.temp_e_s, norm=inpar.norm_s)[0]/de
        flux_h =specx.nthcomp_spec(es=es, gamma=inpar.gamma_h, ktbb=inpar.temp_bb_h, kte=inpar.temp_e_h, norm=inpar.norm_h)[0]/de
        fluxr_d=0.
        fluxr_s=specx.relxillcp_spec(\
            es=es,               incl=inpar.incl,      a=inpar.a,\
            r_in=inpar.r_in_sr,  r_out=inpar.r_out_sr, index=inpar.index_sr,\
            gamma=inpar.gamma_s, logxi=inpar.logxi_sr, logcn=inpar.logcn_sr,\
            cafe=inpar.cafe,     kte=inpar.temp_e_s,   cfref=inpar.cf_sr)[0]/de
        fluxr_h=specx.relxillcp_spec(\
            es=es,               incl=inpar.incl,      a=inpar.a,\
            r_in=inpar.r_in_hr,  r_out=inpar.r_out_hr, index=inpar.index_hr,\
            gamma=inpar.gamma_h, logxi=inpar.logxi_hr, logcn=inpar.logcn_hr,\
            cafe=inpar.cafe,     kte=inpar.temp_e_h,   cfref=inpar.cf_hr)[0]/de
        #print(flux_d, flux_s, flux_h, fluxr_d, fluxr_s, fluxr_h)

        # Normalize such that total corresponds to unity.
        flux_tot=flux_d+flux_s+flux_h+fluxr_d+fluxr_s+fluxr_h
        flux_d /=flux_tot
        flux_s /=flux_tot
        flux_h /=flux_tot
        fluxr_d/=flux_tot
        fluxr_s/=flux_tot
        fluxr_h/=flux_tot
        #print(flux_d, flux_s, flux_h, fluxr_d, fluxr_s, fluxr_h)

        # eta (average)
        # No upper and lower bounds
        eta_d =calc_eta_ave(e_min=inpar.e_minr, e_max=inpar.e_maxr, eta0=inpar.eta0_d, eta1=inpar.eta1_d)
        eta_s =calc_eta_ave(e_min=inpar.e_minr, e_max=inpar.e_maxr, eta0=inpar.eta0_sc, eta1=inpar.eta1_sc)
        eta_h =calc_eta_ave(e_min=inpar.e_minr, e_max=inpar.e_maxr, eta0=inpar.eta0_hc, eta1=inpar.eta1_hc)
        etar_d=0.
        etar_s=calc_eta_ave(e_min=inpar.e_minr, e_max=inpar.e_maxr, eta0=inpar.eta0_sr, eta1=inpar.eta1_sr)
        etar_h=calc_eta_ave(e_min=inpar.e_minr, e_max=inpar.e_maxr, eta0=inpar.eta0_hr, eta1=inpar.eta1_hr)
        #print(eta_d, eta_s, eta_h, etar_d, etar_s, etar_h)

        # \eta*C(E)
        cc_d =flux_d *eta_d
        cc_s =flux_s *eta_s
        cc_h =flux_h *eta_h
        ccr_d=fluxr_d*etar_d
        ccr_s=fluxr_s*etar_s
        ccr_h=fluxr_h*etar_h

        # Energy dissipated
        es_dis=epsilon_calc(r=rs_hcomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
            2.*np.pi*rs_hcomp*wids_hcomp
        # Total energy dissipated (normalization)
        es_dis_tot=np.sum(es_dis)
        # Direct
        ws_hcomp=cc_h*es_dis/es_dis_tot
        # Reprocessed
        wrs_hcomp=ccr_h*es_dis/es_dis_tot

        es_dis=epsilon_calc(r=rs_scomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
            2.*np.pi*rs_scomp*wids_scomp
        es_dis_tot=np.sum(es_dis)
        # Direct
        ws_scomp=cc_s*es_dis/es_dis_tot
        # Reprocessed
        wrs_scomp=ccr_s*es_dis/es_dis_tot

        es_dis=epsilon_calc(r=rs_disk, stress=stress,  gamma=inpar.gamma_disk, r_min=inpar.r_min)*\
            2.*np.pi*rs_disk*wids_disk
        es_dis_tot=np.sum(es_dis)
        # Direct
        ws_disk=cc_d*es_dis/es_dis_tot
        # Reprocessed
        wrs_disk=ccr_d*es_dis/es_dis_tot

        # Direct
        ws_r=ws_disk
        ws_r=np.append(ws_r, ws_scomp)
        ws_r=np.append(ws_r, ws_hcomp)
        # Reprocessed
        wrs_r=wrs_disk
        wrs_r=np.append(wrs_r, wrs_scomp)
        wrs_r=np.append(wrs_r, wrs_hcomp)

        #print(ws_disk)
        #print(wrs_disk)
        #print(ws_r)
        #print(wrs_r)

        # --- Lag --- #
        # Time taken for spectra to respond to mass accretion rate fluctuations
        ### Energy band ###
        # Direct #
        cfs_disk =inpar.cf_disk *np.ones(len(rs_disk))
        cfs_scomp=inpar.cf_scomp*np.ones(len(rs_scomp))
        cfs_hcomp=inpar.cf_hcomp*np.ones(len(rs_hcomp))
        cfs=cfs_disk
        cfs=np.append(cfs, cfs_scomp)
        cfs=np.append(cfs, cfs_hcomp)
        lags=cfs/(rings.fs_vis*bunit.c_rg) # [s]

        # Reprocessed #
        cfs_disk =inpar.cfr_disk *np.ones(len(rs_disk))
        cfs_scomp=inpar.cfr_scomp*np.ones(len(rs_scomp))
        cfs_hcomp=inpar.cfr_hcomp*np.ones(len(rs_hcomp))
        cfs=cfs_disk
        cfs=np.append(cfs, cfs_scomp)
        cfs=np.append(cfs, cfs_hcomp)
        lagrs=cfs/(rings.fs_vis*bunit.c_rg) # [s]

        ### Reference band ###
        # Direct #
        cfs_disk =inpar.cf_diskr *np.ones(len(rs_disk))
        cfs_scomp=inpar.cf_scompr*np.ones(len(rs_scomp))
        cfs_hcomp=inpar.cf_hcompr*np.ones(len(rs_hcomp))
        cfs=cfs_disk
        cfs=np.append(cfs, cfs_scomp)
        cfs=np.append(cfs, cfs_hcomp)
        lags_r=cfs/(rings.fs_vis*bunit.c_rg) # [s]

        # Reprocessed #
        cfs_disk =inpar.cfr_diskr *np.ones(len(rs_disk))
        cfs_scomp=inpar.cfr_scompr*np.ones(len(rs_scomp))
        cfs_hcomp=inpar.cfr_hcompr*np.ones(len(rs_hcomp))
        cfs=cfs_disk
        cfs=np.append(cfs, cfs_scomp)
        cfs=np.append(cfs, cfs_hcomp)
        lagrs_r=cfs/(rings.fs_vis*bunit.c_rg) # [s]

        # --- Impulse response --- #
        # Delay #
        t0s_disk =inpar.t0_disk *np.ones(len(rs_disk))
        t0s_scomp=inpar.t0_scomp*np.ones(len(rs_scomp))
        t0s_hcomp=inpar.t0_hcomp*np.ones(len(rs_hcomp))
        t0s=t0s_disk
        t0s=np.append(t0s, t0s_scomp)
        t0s=np.append(t0s, t0s_hcomp)

        # Duration #
        dt0s_disk =inpar.dt0_disk *np.ones(len(rs_disk))
        dt0s_scomp=inpar.dt0_scomp*np.ones(len(rs_scomp))
        dt0s_hcomp=inpar.dt0_hcomp*np.ones(len(rs_hcomp))
        dt0s=dt0s_disk
        dt0s=np.append(dt0s, dt0s_scomp)
        dt0s=np.append(dt0s, dt0s_hcomp)

        # ----- Print ring information ----- #
        if inpar.display==1:
            print_ring_info(name='R [Rg]',                     xs=rings.rs,                digit=1)
            print_ring_info(name='Disk ring [Rg]',             xs=rs_disk,                 digit=1)
            print_ring_info(name='Soft Compton ring [Rg]',     xs=rs_scomp,                digit=1)
            print_ring_info(name='Hard Compton ring [Rg]',     xs=rs_hcomp,                digit=1)
            print_ring_info(name='Variability frequency [Hz]', xs=rings.fs_vis*bunit.c_rg, digit=3)
            print_ring_info(name='Propagation speed [km/s]',   xs=rings.vs_rad*bunit.c,    digit=3)

        #return

        # ------------------------------------------------------------------------------ #
        # ---------- PSD of mass accretion rate for each ring w/o propagation ---------- #
        # ------------------------------------------------------------------------------ #
        sigs=rings.lfs_var/np.sqrt(rings.n_dec) #\mu is fixed to unity. Probably this does not lose generality. (2021/07/26)
        flupro=FluPro()
        flupro.sigma_set(sigs=sigs)
        flupro.f_set(fs=tifr.fs)
        flupro.f_vis_set(fs_vis=rings.fs_vis)
        flupro.psd_wo_prop(c_rg=bunit.c_rg)

        # ------------------------------------------------------------------------------- #
        # ---------- PSD of mass accretion rate for each ring with propagation ---------- #
        # ------------------------------------------------------------------------------- #
        flupro.dt=tifr.dt
        flupro.n_data=tifr.n_data
        #flupro.psd_w_prop(cds=rings.cds)

        ### 2022/06/04 ###
        ### ------------------------------------------------------ ###
        ### Employ the green function of Rapisarda et al.2017a (3) ###
        ### ------------------------------------------------------ ###
        flupro.psd_w_prop(cs=inpar.cs, fs_prop=rings.fs_acc, dr_r=rings.dr_r, rg_c=bunit.rg_c)

        if inpar.display==1:
            print('--------------------------------------------------')
            print('PSD of the mass accretion rate with propagation was successfully calculated.')

        #return flupro.fs, flupro.psds_prop[-1]

        ##############################################
        ########## Calculate power spectrum ##########
        ##############################################
        if inpar.quant==1:
            # ----------------------- #
            # ----- Energy band ----- #
            # ----------------------- #
            if inpar.display==1:
                print('--------------------------------------------------')
                print('Calculating power spectrum for')
                print('{0:.2f} - {1:.2f} keV'.format(inpar.e_min, inpar.e_max))
            md2fl=Mdot2Flux()
            md2fl.ene_set(e_min=inpar.e_min, e_max=inpar.e_max)
            #################################################################################
            ###### (2021/08/17) Preliminary (haphazard) prescription to set weight, ... #####
            ###### Smarter implementation will be performed.                            #####
            #################################################################################
            # Weight (direct)
            md2fl.ws=ws
            # Weight (reprocessed)
            md2fl.wrs=wrs

            # Lag (direct)
            md2fl.lags=lags
            # Lag (reprocessed)
            md2fl.lagrs=lagrs

            # Impulse response (delay)
            md2fl.t0s=t0s
            # Impulse response (duration)
            md2fl.dt0s=dt0s

            md2fl.mu_fl=1. # <x(E,t)>=\sum _{r_n} w(r_n, E)=1
            md2fl.psd_norm_set(dt=tifr.dt, n_data=tifr.n_data)

            lm2s_prop=flupro.psds_prop/flupro.norm_psd #|M_dot(r, f)|^2

            # fs_vis here is used to calculate accretion time!
            #md2fl.psd_flux_calc(fs=flupro.fs,\
            #                    n_r=rings.n_ring,\
            #                    lm2s=lm2s_prop,\
            #                    fs_vis=rings.fs_acc,\
            #                    cds=rings.cds,\
            #                    xlag=inpar.xlag,\
            #                    dr_r=rings.dr_r,\
            #                    rg_c=bunit.rg_c)
            ### 2022/06/04 ###
            # fs_vis here is used to calculate accretion time!
            md2fl.psd_flux_calc(fs=flupro.fs,\
                                n_r=rings.n_ring,\
                                lm2s=lm2s_prop,\
                                fs_vis=rings.fs_acc,\
                                cs=inpar.cs,\
                                cds=rings.cds,\
                                xlag=inpar.xlag,\
                                dr_r=rings.dr_r,\
                                rg_c=bunit.rg_c)

            psds_fl=md2fl.psd_fl

            if inpar.display==1:
                print('--------------------------------------------------')
                print('PSD of the flux was successfully calculated.')

            ############################################
            ########### Simplest integratinon ##########
            ############################################
            #for i_f_data, f_data in enumerate(fs_data):
            #    i_f=np.abs(flupro.fs-f_data).argmin()
            #    psd=psds_fl[i_f]
            #    if i_f_data==0:
            #        psds=psd
            #    else:
            #        psds=np.append(psds, psd)

            #n_f_data=len(fs_data)
            #for i_f_data in range(n_f_data-1):
            #    # the most simple integration (area of trapezoid)
            #    flux=(psds[i_f_data]+psds[i_f_data+1])*(fs_data[i_f_data+1]-fs_data[i_f_data])/2.
            #    fluxes[i_f_data]=flux

            ################################################
            ########## More accurate integratinon ##########
            ################################################
            n_f_data=len(inpar.fs_data)
            for i_f_data in range(n_f_data-1):
                f_data_min=inpar.fs_data[i_f_data]
                f_data_max=inpar.fs_data[i_f_data+1]
                i_f_min=np.abs(flupro.fs-f_data_min).argmin()
                i_f_max=np.abs(flupro.fs-f_data_max).argmin()
                # In case that frequency bin is too narrow to perform numerical integration
                if i_f_min==i_f_max:
                    psd_int=psds_fl[i_f_min]
                    flux=psd_int*(f_data_max-f_data_min)
                # In case that frequency bin is wide enough to perform numerical integration
                else:
                    fs_int=flupro.fs[i_f_min:i_f_max+1]
                    f_int_min=fs_int[0]
                    f_int_max=fs_int[-1]
                    psds_int=psds_fl[i_f_min:i_f_max+1]
                    # Correction of the difference between the actual integration range and given frequency range
                    flux=integrate.simps(psds_int, fs_int)*(f_data_max-f_data_min)/(f_int_max-f_int_min)
                fluxes[i_f_data]=flux

            #return flupro.fs, psds_fl

        ##############################################
        ########## Calculate cross spectrum ##########
        ##############################################
        elif inpar.quant in [2, 3, 4, 5, 6]:
            if inpar.display==1:
                print('--------------------------------------------------')
                print('Calculating cross spectrum for')
                print('{0:.2f} - {1:.2f} keV vs {2:.2f} - {3:.2f} keV'\
                      .format(inpar.e_min, inpar.e_max, inpar.e_minr, inpar.e_maxr))

            # ----- Weight ----- #
            md2fl=Mdot2Flux()

            ###################
            ### Energy band ###
            ###################
            md2fl.ene_set(e_min=inpar.e_min, e_max=inpar.e_max)
            #################################################################################
            ###### (2021/08/17) Preliminary (haphazard) prescription to set weight, ... #####
            ###### Smarter implementation will be performed.                            #####
            #################################################################################
            # Weight (direct)
            md2fl.ws=ws
            # Weight (reprocessed)
            md2fl.wrs=wrs

            # Lag (direct)
            md2fl.lags=lags
            # Lag (reprocessed)
            md2fl.lagrs=lagrs

            # Impulse response (delay)
            md2fl.t0s=t0s
            # Impulse response (duration)
            md2fl.dt0s=dt0s

            md2fl.mu_fl=1. # <x(E,t)>=\sum _{r_n} w(r_n, E)=1

            ######################
            ### Reference band ###
            ######################
            # Weight (direct)
            md2fl.ws_ref=ws_r
            # Weight (reprocessed)
            md2fl.wrs_ref=wrs_r

            # Lag (direct)
            md2fl.lags_ref=lags_r
            # Lag (reprocessed)
            md2fl.lagrs_ref=lagrs_r

            md2fl.mu_fl_ref=1. # <x(E,t)>=\sum _{r_n} w(r_n, E)=1

            # ----- Normalization ----- #
            md2fl.csd_norm_set(dt=tifr.dt, n_data=tifr.n_data)

            # ----- CSD ----- #
            lm2s_prop=flupro.psds_prop/flupro.norm_psd #|M_dot(r, f)|^2

            # fs_vis here is used to calculate accretion time!
            #md2fl.csd_flux_calc(fs=flupro.fs,\
            #                    n_r=rings.n_ring,\
            #                    lm2s=lm2s_prop,\
            #                    fs_vis=rings.fs_acc,\
            #                    cds=rings.cds,\
            #                    xlag=inpar.xlag,\
            #                    dr_r=rings.dr_r,\
            #                    rg_c=bunit.rg_c)
            ### 2022/06/04 ###
            md2fl.csd_flux_calc(fs=flupro.fs,\
                                n_r=rings.n_ring,\
                                lm2s=lm2s_prop,\
                                fs_vis=rings.fs_acc,\
                                cs=inpar.cs,\
                                cds=rings.cds,\
                                xlag=inpar.xlag,\
                                dr_r=rings.dr_r,\
                                rg_c=bunit.rg_c)

            # Re[CSD]
            if inpar.quant==2:
                csds_fl=np.real(md2fl.csd_fl)
            # Im[CSD]
            elif inpar.quant==3:
                if inpar.invert==1:
                    csds_fl=np.imag(md2fl.csd_fl)
                elif inpar.invert==2:
                    csds_fl=-np.imag(md2fl.csd_fl)
            # |CSD|
            elif inpar.quant==4:
                csds_fl=np.abs(md2fl.csd_fl)
            # Phase lag
            elif inpar.quant==5:
                md2fl.lagf_calc()
                if inpar.invert==1:
                    csds_fl=md2fl.phi
                elif inpar.invert==2:
                    csds_fl=-md2fl.phi
            # Time lag
            elif inpar.quant==6:
                md2fl.lagf_calc()
                if inpar.invert==1:
                    csds_fl=md2fl.tau
                elif inpar.invert==2:
                    csds_fl=-md2fl.tau

            if inpar.display==1:
                print('--------------------------------------------------')
                print('CSD of the flux was successfully calculated.')

            ############################################
            ########### Simplest integratinon ##########
            ############################################
            #for i_f_data, f_data in enumerate(fs_data):
            #    i_f=np.abs(flupro.fs-f_data).argmin()
            #    csd=csds_fl[i_f]
            #    if i_f_data==0:
            #        csds=csd
            #    else:
            #        csds=np.append(csds, csd)

            #n_f_data=len(fs_data)
            #for i_f_data in range(n_f_data-1):
            #    # the most simple integration (area of trapezoid)
            #    flux=(csds[i_f_data]+csds[i_f_data+1])*(fs_data[i_f_data+1]-fs_data[i_f_data])/2.
            #    fluxes[i_f_data]=flux

            ################################################
            ########## More accurate integratinon ##########
            ################################################
            n_f_data=len(inpar.fs_data)
            for i_f_data in range(n_f_data-1):
                f_data_min=inpar.fs_data[i_f_data]
                f_data_max=inpar.fs_data[i_f_data+1]
                i_f_min=np.abs(flupro.fs-f_data_min).argmin()
                i_f_max=np.abs(flupro.fs-f_data_max).argmin()
                # In case that frequency bin is too narrow to perform numerical integration
                if i_f_min==i_f_max:
                    csd_int=csds_fl[i_f_min]
                    flux=csd_int*(f_data_max-f_data_min)
                # In case that frequency bin is wide enough to perform numerical integration
                else:
                    fs_int=flupro.fs[i_f_min:i_f_max+1]
                    f_int_min=fs_int[0]
                    f_int_max=fs_int[-1]
                    csds_int=csds_fl[i_f_min:i_f_max+1]
                    # Correction of the difference between the actual integration range and given frequency range
                    flux=integrate.simps(csds_int, fs_int)*(f_data_max-f_data_min)/(f_int_max-f_int_min)
                fluxes[i_f_data]=flux

            #return flupro.fs, csds_fl

    else:
        print('Error')
        sys.exit()
