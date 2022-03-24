from preflow_h import *
from scipy import integrate

def preflow(engs, params, fluxes):
    # ------------------------------------- #
    # ---------- Initial setting ---------- #
    # ------------------------------------- #
    inpar=SetParameter()
    inpar.set_inpar(pars=params, es=engs)
    #inpar.check_validity()
    
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
    is_hcomp=np.where(rings.rs_min<=inpar.r_mh)[0]
    is_mcomp=np.where((inpar.r_mh<rings.rs_min) & (rings.rs_min<=inpar.r_sm))[0]
    is_scomp=np.where((inpar.r_sm<rings.rs_min) & (rings.rs_min<=inpar.r_ds))[0]
    is_disk =np.where(inpar.r_ds<rings.rs_min)[0]
    rs_flow =rings.rs[is_flow]
    rs_hcomp=rings.rs[is_hcomp]
    rs_mcomp=rings.rs[is_mcomp]
    rs_scomp=rings.rs[is_scomp]
    rs_disk =rings.rs[is_disk]
    # Width is used to calculate weight. (2022/01/28)
    wids_flow =rings.wids[is_flow]
    wids_hcomp=rings.wids[is_hcomp]
    wids_mcomp=rings.wids[is_mcomp]
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

    # --- Damping factor --- #
    cds_flow=(inpar.cd_flow**(1./rings.n_dec))*np.ones(len(rs_flow))
    cds_flow[-1]=inpar.cd_tran
    cds_disk=(inpar.cd_disk**(1./rings.n_dec))*np.ones(len(rs_disk))
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
    rs_mcomp=rs_mcomp[::-1]
    rs_scomp=rs_scomp[::-1]
    rs_disk =rs_disk[::-1]
    wids_flow =wids_flow[::-1]
    wids_hcomp=wids_hcomp[::-1]
    wids_mcomp=wids_mcomp[::-1]
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
    # Energy dissipated
    es_dis=epsilon_calc(r=rs_hcomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
        2.*np.pi*rs_hcomp*wids_hcomp
    # Total energy dissipated (normalization)
    es_dis_tot=np.sum(es_dis)
    # Direct
    ws_hcomp=inpar.cc_hcomp*es_dis/es_dis_tot
    # Reprocessed
    wrs_hcomp=inpar.ccr_hcomp*es_dis/es_dis_tot

    es_dis=epsilon_calc(r=rs_mcomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
        2.*np.pi*rs_mcomp*wids_mcomp
    es_dis_tot=np.sum(es_dis)
    # Direct
    ws_mcomp=inpar.cc_mcomp*es_dis/es_dis_tot
    # Reprocessed
    wrs_mcomp=inpar.ccr_mcomp*es_dis/es_dis_tot

    es_dis=epsilon_calc(r=rs_scomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
        2.*np.pi*rs_scomp*wids_scomp
    es_dis_tot=np.sum(es_dis)
    # Direct
    ws_scomp=inpar.cc_scomp*es_dis/es_dis_tot
    # Reprocessed
    wrs_scomp=inpar.ccr_scomp*es_dis/es_dis_tot

    es_dis=epsilon_calc(r=rs_disk, stress=stress,  gamma=inpar.gamma_disk, r_min=inpar.r_min)*\
        2.*np.pi*rs_disk*wids_disk
    es_dis_tot=np.sum(es_dis)
    # Direct
    ws_disk=inpar.cc_disk*es_dis/es_dis_tot
    # Reprocessed
    wrs_disk=inpar.ccr_disk*es_dis/es_dis_tot

    # Direct
    ws=ws_disk
    ws=np.append(ws, ws_scomp)
    ws=np.append(ws, ws_mcomp)
    ws=np.append(ws, ws_hcomp)
    # Reprocessed
    wrs=wrs_disk
    wrs=np.append(wrs, wrs_scomp)
    wrs=np.append(wrs, wrs_mcomp)
    wrs=np.append(wrs, wrs_hcomp)


    ### Reference band ###
    # Energy dissipated
    es_dis=epsilon_calc(r=rs_hcomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
        2.*np.pi*rs_hcomp*wids_hcomp
    # Total energy dissipated (normalization)
    es_dis_tot=np.sum(es_dis)
    # Direct
    ws_hcomp=inpar.cc_hcompr*es_dis/es_dis_tot
    # Reprocessed
    wrs_hcomp=inpar.ccr_hcompr*es_dis/es_dis_tot

    es_dis=epsilon_calc(r=rs_mcomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
        2.*np.pi*rs_mcomp*wids_mcomp
    es_dis_tot=np.sum(es_dis)
    # Direct
    ws_mcomp=inpar.cc_mcompr*es_dis/es_dis_tot
    # Reprocessed
    wrs_mcomp=inpar.ccr_mcompr*es_dis/es_dis_tot

    es_dis=epsilon_calc(r=rs_scomp, stress=stress, gamma=inpar.gamma_flow, r_min=inpar.r_min)*\
        2.*np.pi*rs_scomp*wids_scomp
    es_dis_tot=np.sum(es_dis)
    # Direct
    ws_scomp=inpar.cc_scompr*es_dis/es_dis_tot
    # Reprocessed
    wrs_scomp=inpar.ccr_scompr*es_dis/es_dis_tot

    es_dis=epsilon_calc(r=rs_disk, stress=stress,  gamma=inpar.gamma_disk, r_min=inpar.r_min)*\
        2.*np.pi*rs_disk*wids_disk
    es_dis_tot=np.sum(es_dis)
    # Direct
    ws_disk=inpar.cc_diskr*es_dis/es_dis_tot
    # Reprocessed
    wrs_disk=inpar.ccr_diskr*es_dis/es_dis_tot

    # Direct
    ws_r=ws_disk
    ws_r=np.append(ws_r, ws_scomp)
    ws_r=np.append(ws_r, ws_mcomp)
    ws_r=np.append(ws_r, ws_hcomp)
    # Reprocessed
    wrs_r=wrs_disk
    wrs_r=np.append(wrs_r, wrs_scomp)
    wrs_r=np.append(wrs_r, wrs_mcomp)
    wrs_r=np.append(wrs_r, wrs_hcomp)
    
    # --- Lag --- #
    # Time taken for spectra to respond to mass accretion rate fluctuations
    ### Energy band ###
    # Direct #
    cfs_disk =inpar.cf_disk *np.ones(len(rs_disk))
    cfs_scomp=inpar.cf_scomp*np.ones(len(rs_scomp))
    cfs_mcomp=inpar.cf_mcomp*np.ones(len(rs_mcomp))
    cfs_hcomp=inpar.cf_hcomp*np.ones(len(rs_hcomp))
    cfs=cfs_disk
    cfs=np.append(cfs, cfs_scomp)
    cfs=np.append(cfs, cfs_mcomp)
    cfs=np.append(cfs, cfs_hcomp)
    lags=cfs/(rings.fs_vis*bunit.c_rg) # [s]

    # Reprocessed #
    cfs_disk =inpar.cfr_disk *np.ones(len(rs_disk))
    cfs_scomp=inpar.cfr_scomp*np.ones(len(rs_scomp))
    cfs_mcomp=inpar.cfr_mcomp*np.ones(len(rs_mcomp))
    cfs_hcomp=inpar.cfr_hcomp*np.ones(len(rs_hcomp))
    cfs=cfs_disk
    cfs=np.append(cfs, cfs_scomp)
    cfs=np.append(cfs, cfs_mcomp)
    cfs=np.append(cfs, cfs_hcomp)
    lagrs=cfs/(rings.fs_vis*bunit.c_rg) # [s]

    ### Reference band ###
    # Direct #
    cfs_disk =inpar.cf_diskr *np.ones(len(rs_disk))
    cfs_scomp=inpar.cf_scompr*np.ones(len(rs_scomp))
    cfs_mcomp=inpar.cf_mcompr*np.ones(len(rs_mcomp))
    cfs_hcomp=inpar.cf_hcompr*np.ones(len(rs_hcomp))
    cfs=cfs_disk
    cfs=np.append(cfs, cfs_scomp)
    cfs=np.append(cfs, cfs_mcomp)
    cfs=np.append(cfs, cfs_hcomp)
    lags_r=cfs/(rings.fs_vis*bunit.c_rg) # [s]

    # Reprocessed #
    cfs_disk =inpar.cfr_diskr *np.ones(len(rs_disk))
    cfs_scomp=inpar.cfr_scompr*np.ones(len(rs_scomp))
    cfs_mcomp=inpar.cfr_mcompr*np.ones(len(rs_mcomp))
    cfs_hcomp=inpar.cfr_hcompr*np.ones(len(rs_hcomp))
    cfs=cfs_disk
    cfs=np.append(cfs, cfs_scomp)
    cfs=np.append(cfs, cfs_mcomp)
    cfs=np.append(cfs, cfs_hcomp)
    lagrs_r=cfs/(rings.fs_vis*bunit.c_rg) # [s]

    # --- Impulse response --- #
    # Delay #
    t0s_disk =inpar.t0_disk *np.ones(len(rs_disk))
    t0s_scomp=inpar.t0_scomp*np.ones(len(rs_scomp))
    t0s_mcomp=inpar.t0_mcomp*np.ones(len(rs_mcomp))
    t0s_hcomp=inpar.t0_hcomp*np.ones(len(rs_hcomp))
    t0s=t0s_disk
    t0s=np.append(t0s, t0s_scomp)
    t0s=np.append(t0s, t0s_mcomp)
    t0s=np.append(t0s, t0s_hcomp)

    # Duration #
    dt0s_disk =inpar.dt0_disk *np.ones(len(rs_disk))
    dt0s_scomp=inpar.dt0_scomp*np.ones(len(rs_scomp))
    dt0s_mcomp=inpar.dt0_mcomp*np.ones(len(rs_mcomp))
    dt0s_hcomp=inpar.dt0_hcomp*np.ones(len(rs_hcomp))
    dt0s=dt0s_disk
    dt0s=np.append(dt0s, dt0s_scomp)
    dt0s=np.append(dt0s, dt0s_mcomp)
    dt0s=np.append(dt0s, dt0s_hcomp)

    # ----- Print ring information ----- #
    if inpar.display==1:
        print_ring_info(name='R [Rg]',                     xs=rings.rs,                digit=1)
        print_ring_info(name='Disk ring [Rg]',             xs=rs_disk,                 digit=1)
        print_ring_info(name='Soft Compton ring [Rg]',     xs=rs_scomp,                digit=1)
        print_ring_info(name='Mid Compton ring [Rg]',      xs=rs_mcomp,                digit=1)
        print_ring_info(name='Hard Compton ring [Rg]',     xs=rs_hcomp,                digit=1)
        print_ring_info(name='Variability frequency [Hz]', xs=rings.fs_vis*bunit.c_rg, digit=3)
        #print_ring_info(name='Propagation frequency [Hz]', xs=rings.fs_acc*bunit.c_rg, digit=3)
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
    flupro.psd_w_prop(cds=rings.cds)
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
        md2fl.psd_flux_calc(fs=flupro.fs,\
                            n_r=rings.n_ring,\
                            lm2s=lm2s_prop,\
                            fs_vis=rings.fs_acc,\
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
        md2fl.csd_flux_calc(fs=flupro.fs,\
                            n_r=rings.n_ring,\
                            lm2s=lm2s_prop,\
                            fs_vis=rings.fs_acc,\
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
