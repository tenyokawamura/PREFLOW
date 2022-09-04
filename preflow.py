from preflow_h import *
from scipy import integrate

def preflow(engs, params, fluxes):
    # ------------------------------------- #
    # ---------- Initial setting ---------- #
    # ------------------------------------- #
    inpar=SetParameter()
    inpar.preflow_set_inpar(pars=params, es=engs)
    
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
    rings.r_sh =inpar.r_sh
    rings.r_ds =inpar.r_ds
    rings.r_out=inpar.r_out
    rings.t0_d =inpar.t0_d
    rings.t0_s =inpar.t0_s
    rings.t0_h =inpar.t0_h
    rings.dt0_d=inpar.dt0_d
    rings.dt0_s=inpar.dt0_s
    rings.dt0_h=inpar.dt0_h
    rings.par_set()
    rings.coeff_calc()
    rings.lowboundr_calc()
    rings.centerr_calc()
    rings.width_calc()
    rings.interval_calc()
    rings.n_ring_dec_calc()
    rings.ring_classify()
    rings.variability_frequency_calc(
        cb_f=inpar.cb_f, m_f=inpar.m_f, cbp_f=inpar.cbp_f, mp_f=inpar.mp_f,\
        cb_d=inpar.cb_d, m_d=inpar.m_d, cbp_d=inpar.cbp_d, mp_d=inpar.mp_d)
    rings.damping_calc()
    rings.spec_lag_calc()
    rings.variability_amplitude_calc(\
        cf_var_f=inpar.cf_var_f, dr_var_f=inpar.dr_var_f,\
        cf_var_d=inpar.cf_var_d, dr_var_d=inpar.dr_var_d)
    rings.imp_resp_set()

    # --------------------------------------- #
    # ---------- Calculate weight  ---------- #
    # --------------------------------------- #
    # --- No reverberation --- #
    if [inpar.cs_sr, inpar.cs_hr, inpar.csp_sr, inpar.csp_hr,\
        inpar.cs_sr_r, inpar.cs_hr_r, inpar.csp_sr_r, inpar.csp_hr_r]==\
       [0., 0., 0., 0., 0., 0., 0., 0.]:
        # Subject band
        rings.ws=preflow_weight_calc(\
            cs_d=inpar.cs_d, cs_s=inpar.cs_s, cs_h=inpar.cs_h,\
            csp_s=inpar.csp_s, csp_h=inpar.csp_h,\
            rs_d=rings.rs_d,       rs_s=rings.rs_s,       rs_h=rings.rs_h,\
            wids_d=rings.wids_d,   wids_s=rings.wids_s,   wids_h=rings.wids_h,\
            stress=inpar.stress,   r_min=inpar.r_min,     index_d=inpar.gamma_disk, index_f=inpar.gamma_flow)
        # Reference band
        rings.ws_r=preflow_weight_calc(\
            cs_d=inpar.cs_d_r, cs_s=inpar.cs_s_r, cs_h=inpar.cs_h_r,\
            csp_s=inpar.csp_s_r, csp_h=inpar.csp_h_r,\
            rs_d=rings.rs_d,       rs_s=rings.rs_s,       rs_h=rings.rs_h,\
            wids_d=rings.wids_d,   wids_s=rings.wids_s,   wids_h=rings.wids_h,\
            stress=inpar.stress,   r_min=inpar.r_min,     index_d=inpar.gamma_disk, index_f=inpar.gamma_flow)
        rings.wrss=[0.]   # Unused
        rings.wrss_r=[0.] # Unused

    # --- Reverberation taken into account --- #
    else:
        # Subject band
        rings.ws, rings.wrss=preflow_ref_weight_calc(\
            cs_d=inpar.cs_d, cs_s=inpar.cs_s, cs_h=inpar.cs_h, cs_sr=inpar.cs_sr, cs_hr=inpar.cs_hr,\
            csp_s=inpar.csp_s, csp_h=inpar.csp_h, csp_sr=inpar.csp_sr, csp_hr=inpar.csp_hr,\
            rs_d=rings.rs_d,         rs_s=rings.rs_s,       rs_h=rings.rs_h,\
            wids_d=rings.wids_d,     wids_s=rings.wids_s,   wids_h=rings.wids_h,\
            stress=inpar.stress,     r_min=inpar.r_min,     index_d=inpar.gamma_disk, index_f=inpar.gamma_flow,\
            t0_d=inpar.t0_d,         dt0_d=inpar.dt0_d,\
            t0_s=inpar.t0_s,         dt0_s=inpar.dt0_s,\
            t0_h=inpar.t0_h,         dt0_h=inpar.dt0_h,\
            fs=tifr.fs)
        # Reference band
        rings.ws_r, rings.wrss_r=preflow_ref_weight_calc(\
            cs_d=inpar.cs_d_r, cs_s=inpar.cs_s_r, cs_h=inpar.cs_h_r, cs_sr=inpar.cs_sr_r, cs_hr=inpar.cs_hr_r,\
            csp_s=inpar.csp_s_r, csp_h=inpar.csp_h_r, csp_sr=inpar.csp_sr_r, csp_hr=inpar.csp_hr_r,\
            rs_d=rings.rs_d,         rs_s=rings.rs_s,       rs_h=rings.rs_h,\
            wids_d=rings.wids_d,     wids_s=rings.wids_s,   wids_h=rings.wids_h,\
            stress=inpar.stress,     r_min=inpar.r_min,     index_d=inpar.gamma_disk, index_f=inpar.gamma_flow,\
            t0_d=inpar.t0_d,         dt0_d=inpar.dt0_d,\
            t0_s=inpar.t0_s,         dt0_s=inpar.dt0_s,\
            t0_h=inpar.t0_h,         dt0_h=inpar.dt0_h,\
            fs=tifr.fs)

    # --- Final operation --- #
    rings.out2in()

    # ----- Print ring information ----- #
    if inpar.display==1:
        print_ring_info(name='R [Rg]',                     xs=rings.rs,                digit=1)
        print_ring_info(name='Disk ring [Rg]',             xs=rings.rs_d,              digit=1)
        print_ring_info(name='Soft Compton ring [Rg]',     xs=rings.rs_s,              digit=1)
        print_ring_info(name='Hard Compton ring [Rg]',     xs=rings.rs_h,              digit=1)
        print_ring_info(name='Variability frequency [Hz]', xs=rings.fs_vis*bunit.c_rg, digit=3)
        print_ring_info(name='Propagation speed [km/s]',   xs=rings.vs_prop*bunit.c,   digit=3)

    #print(rings.ws)
    #print(rings.wrs)
    #print(rings.ws_r)
    #print(rings.wrs_r)
    #return

    # ----------------------------------------------------------- #
    # ---------- Power spectrum of mass accretion rate ---------- #
    # ----------------------------------------------------------- #
    sigs=rings.lfs_var/np.sqrt(rings.n_dec) #\mu is fixed to unity. Probably this does not lose generality. (2021/07/26)
    flupro=FluPro()
    flupro.sigma_set(sigs=sigs)
    flupro.f_set(fs=tifr.fs)
    flupro.f_vis_set(fs_vis=rings.fs_vis)
    flupro.psd_wo_prop(c_rg=bunit.c_rg)

    flupro.dt=tifr.dt
    flupro.n_data=tifr.n_data
    #flupro.psd_w_prop() # No damping
    ### 2022/06/04 ###
    ### Employ the green function of Rapisarda et al.2017a (3) ###
    flupro.psd_w_prop(cs=inpar.cs, fs_prop=rings.fs_prop, dr_r=rings.dr_r, rg_c=bunit.rg_c) # With damping

    if inpar.display==1:
        print('--------------------------------------------------')
        print('PSD of the mass accretion rate with propagation was successfully calculated.')

    #return flupro.fs, flupro.psds_prop[-1]
    #return

    # -------------------------------------------------- #
    # ---------- Power/Cross spectrum of flux ---------- #
    # -------------------------------------------------- #
    md2fl=Mdot2Flux()
    md2fl.ene_set(e_min=inpar.e_min, e_max=inpar.e_max)
    # Number of rings
    md2fl.n_ring=rings.n_ring
    # Propagation frequency
    md2fl.fs_prop=rings.fs_prop
    # Damping (Smoothing) factor
    md2fl.cs=inpar.cs
    # dr/r
    md2fl.dr_r=rings.dr_r
    # Rg/c
    md2fl.rg_c=bunit.rg_c
    # Impulse response (delay)
    md2fl.t0s=rings.t0s
    # Impulse response (duration)
    md2fl.dt0s=rings.dt0s
    # Timing interval
    md2fl.dt=tifr.dt
    # Number of data
    md2fl.n_data=tifr.n_data
    # Power spectrum of mass accretion rate
    md2fl.fs=flupro.fs
    md2fl.lm2s=flupro.psds_prop/flupro.norm_psd #|M_dot(r, f)|^2

    del flupro

    # ----- Power spectrum ----- #
    if inpar.quant==1:
        # ----------------------- #
        # ----- Energy band ----- #
        # ----------------------- #
        if inpar.display==1:
            print('--------------------------------------------------')
            print('Calculating power spectrum for')
            print('{0:.2f} - {1:.2f} keV'.format(inpar.e_min, inpar.e_max))

        # --- No reverberation --- #
        if [inpar.cs_sr, inpar.cs_hr, inpar.csp_sr, inpar.csp_hr,\
            inpar.cs_sr_r, inpar.cs_hr_r, inpar.csp_sr_r, inpar.csp_hr_r]==\
           [0., 0., 0., 0., 0., 0., 0., 0.]:
            # Weight (direct)
            md2fl.ws=rings.ws
            # Time-averaged flux
            md2fl.mu_fl=1. # <x(E,t)>=1

            md2fl.psd_norm_set()
            md2fl.psd_flux_calc()

        # --- Reverberation taken into account --- #
        else:
            # Weight (direct)
            md2fl.ws=rings.ws
            # Weight (reprocessed)
            md2fl.wrss=rings.wrss
            # Time-averaged flux
            md2fl.mu_fl=1. # <x(E,t)>=1

            md2fl.psd_norm_set()
            md2fl.psd_flux_ref_calc()

        #psds_fl=md2fl.psd_fl
        data=md2fl.psd_fl

        if inpar.display==1:
            print('--------------------------------------------------')
            print('PSD of the flux was successfully calculated.')

        #return md2fl.fs, md2fl.psds_fl

    # ----- Cross spectrum of flux ----- #
    elif inpar.quant in [2, 3, 4, 5, 6]:
        if inpar.display==1:
            print('--------------------------------------------------')
            print('Calculating cross spectrum for')
            print('{0:.2f} - {1:.2f} keV vs {2:.2f} - {3:.2f} keV'\
                  .format(inpar.e_min, inpar.e_max, inpar.e_minr, inpar.e_maxr))

        # --- No reverberation --- #
        if [inpar.cs_sr, inpar.cs_hr, inpar.csp_sr, inpar.csp_hr,\
            inpar.cs_sr_r, inpar.cs_hr_r, inpar.csp_sr_r, inpar.csp_hr_r]==\
           [0., 0., 0., 0., 0., 0., 0., 0.]:
            # --- Subject band --- #
            # Weight (direct)
            md2fl.ws=rings.ws
            # Time-averaged flux
            md2fl.mu_fl=1. # <x(E,t)>=1

            # --- Reference band --- #
            # Weight (direct)
            md2fl.ws_r=rings.ws_r
            # Time-averaged flux
            md2fl.mu_fl_r=1. # <x(E,t)>=1

            md2fl.csd_norm_set()
            md2fl.csd_flux_calc()

        # --- Reverberation taken into account --- #
        else:
            # --- Subject band --- #
            # Weight (direct)
            md2fl.ws=rings.ws
            # Weight (reprocessed)
            md2fl.wrss=rings.wrss
            # Time-averaged flux
            md2fl.mu_fl=1. # <x(E,t)>=1

            # --- Reference band --- #
            # Weight (direct)
            md2fl.ws_r=rings.ws_r
            # Weight (reprocessed)
            md2fl.wrss_r=rings.wrss_r
            # Time-averaged flux
            md2fl.mu_fl_r=1. # <x(E,t)>=1

            md2fl.csd_norm_set()
            md2fl.csd_flux_ref_calc()

        # Choice of quantity
        md2fl.quant=inpar.quant
        # Choice of polarity
        md2fl.invert=inpar.invert

        md2fl.csd_convert()

        #csds_fl=md2fl.csd_fl
        data=md2fl.csd_fl

        if inpar.display==1:
            print('--------------------------------------------------')
            print('CSD of the flux was successfully calculated.')

        #return md2fl.fs, md2fl.csds_fl

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
        i_f_min=np.abs(md2fl.fs-f_data_min).argmin()
        i_f_max=np.abs(md2fl.fs-f_data_max).argmin()
        # In case that frequency bin is too narrow to perform numerical integration
        if i_f_min==i_f_max:
            datum_int=data[i_f_min]
            flux=datum_int*(f_data_max-f_data_min)
        # In case that frequency bin is wide enough to perform numerical integration
        else:
            fs_int=md2fl.fs[i_f_min:i_f_max+1]
            f_int_min=fs_int[0]
            f_int_max=fs_int[-1]
            data_int=data[i_f_min:i_f_max+1]
            # Correction of the difference between the actual integration range and given frequency range
            flux=integrate.simps(data_int, fs_int)*(f_data_max-f_data_min)/(f_int_max-f_int_min)
        fluxes[i_f_data]=flux
