from preflow_h import *
from scipy import integrate

def preflow(engs, params, fluxes):
    # ------------------------------------- #
    # ---------- Initial setting ---------- #
    # ------------------------------------- #
    mass        =params[0]  # BH mass [solar mass]
    r_in        =params[1]  # Inner radius of hard Compton (inner radius of hot flow) [Rg]
    r_sh        =params[2]  # Transition radius between hard and soft Compton [Rg]
    r_ds        =params[3]  # Transition radius between soft Compton and variable disk [Rg]
    r_out       =params[4]  # Transition radius between soft Compton and variable disk (outer radius of hot flow) [Rg]
    n_ring      =params[5]  # Outer radius of variable disk [Rg]
    tref        =params[6]  # Start time of reflection impulse response [sec]
    dtref       =params[7]  # Time width of reflection impulse response [sec]
    lf_var      =params[8]  # Fractional variability of mass accretion rate in radial decade [-]
    lb_disk     =params[9]  # B_{disk} [-]
    m_disk      =params[10] # m_{disk} [-]
    lb_flow     =params[11] # B_{flow} [-]
    m_flow      =params[12] # m_{flow} [-]
    stress      =params[13] # 1: stressed, 2: stress-free in emissivity
    gamma       =params[14] # Radial index of emissivity [-]
    e_min       =params[15] # Lower bound of energy band [keV] (unused)
    e_max       =params[16] # Upper bound of energy band [keV] (unused)
    frac_disk   =params[17] # Fraction of variable disk in the energy band [counts keV^-1 s^-1]
    frac_scomp  =params[18] # Fraction of soft Compton in the energy band [counts keV^-1 s^-1]
    frac_hcomp  =params[19] # Fraction of hard Compton in the energy band [counts keV^-1 s^-1]
    frac_sref   =params[20] # Fraction of soft reflection in the energy band [counts keV^-1 s^-1]
    frac_href   =params[21] # Fraction of hard reflection in the energy band [counts keV^-1 s^-1]
    e_minr      =params[22] # Lower bound of reference band [keV] (unused)
    e_maxr      =params[23] # Upper bound of reference band [keV] (unused)
    frac_diskr  =params[24] # Fraction of variable disk in the reference band [counts keV^-1 s^-1]
    frac_scompr =params[25] # Fraction of soft Compton in the reference band [counts keV^-1 s^-1]
    frac_hcompr =params[26] # Fraction of hard Compton in the reference band [counts keV^-1 s^-1]
    frac_srefr  =params[27] # Fraction of soft reflection in the reference band [counts keV^-1 s^-1]
    frac_hrefr  =params[28] # Fraction of hard reflection in the reference band [counts keV^-1 s^-1]
    e_minrr     =params[29] # Lower bound of reference band 'for reflection' [keV] (unused)
    e_maxrr     =params[30] # Upper bound of reference band 'for reflection' [keV] (unused)
    frac_scomprr=params[31] # Soft Compton in the reference band 'for reflection' [counts kev^-1 s^-1]
    quant       =params[32]
        # 1: power spectrum 
        # 2: real part of cross spectrum
        # 3: imaginary part of cross spectrum
        # 4: absolute value of cross spectrum
        # 5: phase lag (Positive lag means reference band logging behind energy band.)
        # 6: time lag  (Positive lag means reference band logging behind energy band.)
    display     =params[33] # 1: display, 2: not display
    
    # Parameters, which are no longer free.
    frac_tot =1. # Total fraction in the energy band [counts keV^-1 s^-1], i.e., 1.
    frac_totr=1. # Total fraction in the reference band [counts keV^-1 s^-1], 1.e., 1.
    frac_hcomprr=1.-frac_scomprr # Fraction of hard Compton in the hot flow in the reference band 'for reflection' [counts kev^-1 s^-1]

    # Impulse response of reflection
    t0=tref+(dtref/2.)
    dt0=dtref

    # PREFLOW model is a timing model!
    # Energy in XSPEC corresponds to Fourier frequency in preflow.
    fs_data=engs

    # --------------------------------------- #
    # ---------- Check constraints ---------- #
    # --------------------------------------- #
    ############
    ### Flux ###
    ############
    # Energy band
    if frac_disk+frac_scomp+frac_hcomp+frac_sref+frac_href<=frac_tot:
        pass
    elif frac_disk>frac_tot:
        print('Error: frac_disk>1 --> frac_disk=1, frac_scomp=frac_hcomp=frac_sref=frac_href=0')
        frac_disk=frac_tot
        frac_scomp=0
        frac_hcomp=0
        frac_sref=0
        frac_href=0
    elif frac_disk+frac_scomp>frac_tot:
        print('Error: frac_disk+frac_scomp>1 --> frac_scomp=1-frac_disk, frac_hcomp=frac_sref=frac_href=0')
        frac_scomp=frac_tot-frac_disk
        frac_hcomp=0
        frac_sref=0
        frac_href=0
    elif frac_disk+frac_scomp+frac_hcomp>frac_tot:
        print('Error: frac_disk+frac_scomp+frac_hcomp>1 --> frac_hcomp=1-(frac_disk+frac_scomp), frac_sref=frac_href=0')
        frac_hcomp=frac_tot-(frac_disk+frac_scomp)
        frac_sref=0
        frac_href=0
    elif frac_disk+frac_scomp+frac_hcomp+frac_sref>frac_tot:
        print('Error: frac_disk+frac_scomp+frac_hcomp+frac_sref>1 --> frac_sref=1-(frac_disk+frac_scomp+frac_hcomp), frac_href=0')
        frac_sref=frac_tot-(frac_disk+frac_scomp+frac_hcomp)
        frac_href=0
    elif frac_disk+frac_scomp+frac_hcomp+frac_sref+frac_href>frac_tot:
        print('Error: frac_disk+frac_scomp+frac_hcomp+frac_sref>1 --> frac_href=1-(frac_disk+frac_scomp+frac_hcomp+frac_sref)')
        frac_href=frac_tot-(frac_disk+frac_scomp+frac_hcomp+frac_sref)
    else:
        print('Error')

    # Reference band
    if frac_diskr+frac_scompr+frac_hcompr+frac_srefr+frac_hrefr<=frac_tot:
        pass
    elif frac_diskr>frac_tot:
        print('Error: frac_diskr>1 --> frac_diskr=1, frac_scompr=frac_hcompr=frac_srefr=frac_hrefr=0')
        frac_diskr=frac_tot
        frac_scompr=0.
        frac_hcompr=0.
        frac_srefr=0.
        frac_hrefr=0.
    elif frac_diskr+frac_scompr>frac_tot:
        print('Error: frac_diskr+frac_scompr>1 --> frac_scompr=1-frac_diskr, frac_hcompr=frac_srefr=frac_hrefr=0')
        frac_scompr=frac_tot-frac_diskr
        frac_hcompr=0.
        frac_srefr=0.
        frac_hrefr=0.
    elif frac_diskr+frac_scompr+frac_hcompr>frac_tot:
        print('Error: frac_diskr+frac_scompr+frac_hcompr>1 --> frac_hcompr=1-(frac_diskr+frac_scompr), frac_srefr=frac_hrefr=0')
        frac_hcompr=frac_tot-(frac_diskr+frac_scompr)
        frac_srefr=0.
        frac_hrefr=0.
    elif frac_diskr+frac_scompr+frac_hcompr+frac_srefr>frac_tot:
        print('Error: frac_diskr+frac_scompr+frac_hcompr+frac_srefr>1 --> frac_srefr=1-(frac_diskr+frac_scompr+frac_hcompr), frac_hrefr=0')
        frac_srefr=frac_tot-(frac_diskr+frac_scompr+frac_hcompr)
        frac_hrefr=0.
    elif frac_diskr+frac_scompr+frac_hcompr+frac_srefr+frac_hrefr>frac_tot:
        print('Error: frac_diskr+frac_scompr+frac_hcompr+frac_srefr>1 --> frac_hrefr=1-(frac_diskr+frac_scompr+frac_hcompr+frac_srefr)')
        frac_hrefr=frac_tot-(frac_diskr+frac_scompr+frac_hcompr+frac_srefr)
    else:
        print('Error')

    # Reference band 'for reflection'
    if frac_scomprr<=1.:
        pass
    else:
        print('Error: frac_scomprr>1 --> frac_scomprr=1, frac_hcomprr=0')
        frac_scomprr=1.
        frac_hcomprr=0.

    ################
    ### Geometry ###
    ################
    if r_in>r_sh:
        print('Error: rin>rsh --> rsh=rin')
        r_sh=r_in
    if r_sh>r_ds:
        print('Error: rsh>rds --> rds=rsh')
        r_ds=sh
    if r_ds>r_out:
        print('Error: rds>rout --> rout=rds')
        r_ds=r_out

    # ------------------------------------ #
    # ---------- Set basic unit ---------- #
    # ------------------------------------ #
    bunit=BasicUnit()
    bunit.set_unit(mass=mass)
    # ----- Print basic information ----- #
    if display==1:
        print('--------------------------------------------------')
        print('M   : {0:.1f} solar mass'.format(bunit.m))
        print('Rg  : {0:.1f} km'.format(bunit.rg))
        #print('Rg/c: {0:.1e} sec'.format(bunit.rg_c))
        #print('c/Rg: {0:.1e} Hz'.format(bunit.c_rg))

    # ---------------------------------------- #
    # ---------- Set time/frequency ---------- #
    # ---------------------------------------- #
    f_data_min=fs_data[0]
    f_data_max=fs_data[-1]

    # ----- Decide dt ----- #
    # dt = coeff x 10^{index}
    coeff=2
    index=1
    while True:
        if coeff==1:
            coeff=5
            index-=1
        elif coeff==2:
            coeff=1
        elif coeff==5:
            coeff=2
        else:
            print('Error')
            sys.exit()
        dt=coeff*10**(index)
        f_max=1./(2.*dt)
        if f_data_max<f_max:
            break

    # ----- Decide n_data ----- #
    n_data=1
    while True:
        n_data*=2
        f_min=1./(n_data*dt)
        if f_min<f_data_min:
            break
    
    tifr=TimeFrequency()
    tifr.dt=dt
    tifr.n_data=n_data
    tifr.t_set()
    tifr.f_set()

    # ----- Print timing information ----- #
    if display==1:
        print('--------------------------------------------------')
        print('Timing interval: {0:.0e} sec'.format(tifr.dt))
        print('Number of data: {0:.0f}'.format(tifr.n_data))

    # ----------------------------------------------- #
    # ---------- Accretion flow separation ---------- #
    # ----------------------------------------------- #
    rings=Flow2Ring()
    rings.n_ring=int(n_ring)
    rings.r_in=r_in
    rings.r_out=r_out
    rings.par_set()
    rings.coeff_calc()
    rings.lowboundr_calc()
    rings.centerr_calc()
    rings.width_calc()
    rings.interval_calc()
    rings.n_ring_dec_calc()

    # ------------------------------- #
    # ---------- Set flux  ---------- #
    # ------------------------------- #
    spec=FluxData()
    spec.set_flux(e_min  =e_min,\
                  e_max  =e_max,\
                  disk   =frac_disk,\
                  scomp  =frac_scomp,\
                  hcomp  =frac_hcomp,\
                  sref   =frac_sref,\
                  href   =frac_href,\
                  tot    =frac_tot,\
                  e_minr =e_minr,\
                  e_maxr =e_maxr,\
                  diskr  =frac_diskr,\
                  scompr =frac_scompr,\
                  hcompr =frac_hcompr,\
                  srefr  =frac_srefr,\
                  hrefr  =frac_hrefr,\
                  totr   =frac_totr,\
                  e_minrr=e_minrr,\
                  e_maxrr=e_maxrr,\
                  scomprr=frac_scomprr,\
                  hcomprr=frac_hcomprr)

    # ------------------------------------------ #
    # ---------- Assign flux to ring  ---------- #
    # ------------------------------------------ #
    r_min_eps=r_in
    ###########################
    ### Hard Comptonization ###
    ###########################
    hcomp=Ring2Spec()
    hcomp.r_min=r_in
    hcomp.r_max=r_sh
    hcomp.lb=lb_flow
    hcomp.m=m_flow
    hcomp.epsilon_par_set(stress=stress, gamma=gamma, r_min=r_min_eps)
    hcomp.ring_assign(rs=rings.rs,\
                      rs_min=rings.rs_min,\
                      wids=rings.wids,\
                      drs=rings.drs)
    hcomp.f_vis_set()
    hcomp.v_rad_set()
    hcomp.epsilon_set()

    ###########################
    ### Soft Comptonization ###
    ###########################
    scomp=Ring2Spec()
    scomp.r_min=r_sh
    scomp.r_max=r_ds
    scomp.lb=lb_flow
    scomp.m=m_flow
    scomp.epsilon_par_set(stress=stress, gamma=gamma, r_min=r_min_eps)
    scomp.ring_assign(rs=rings.rs,\
                      rs_min=rings.rs_min,\
                      wids=rings.wids,\
                      drs=rings.drs)
    scomp.f_vis_set()
    scomp.v_rad_set()
    scomp.epsilon_set()

    ######################
    ### Disk blackbody ###
    ######################
    disk=Ring2Spec()
    disk.r_min=r_ds
    disk.r_max=r_out
    disk.lb=lb_disk
    disk.m=m_disk
    disk.epsilon_par_set(stress=stress, gamma=gamma, r_min=r_min_eps)
    disk.ring_assign(rs=rings.rs,\
                     rs_min=rings.rs_min,\
                     wids=rings.wids,\
                     drs=rings.drs)
    disk.f_vis_set()
    disk.v_rad_set()
    disk.epsilon_set()

    # --- Give info on viscous frequency and radial velocity to Flow2Ring class instance --- #
    fs_vis=hcomp.fs_vis
    fs_vis=np.append(fs_vis, scomp.fs_vis)
    fs_vis=np.append(fs_vis, disk.fs_vis)
    rings.fs_vis=fs_vis

    vs_rad=hcomp.vs_rad
    vs_rad=np.append(vs_rad, scomp.vs_rad)
    vs_rad=np.append(vs_rad, disk.vs_rad)
    rings.vs_rad=vs_rad

    eps=hcomp.eps
    eps=np.append(eps, scomp.eps)
    eps=np.append(eps, disk.eps)
    rings.eps=eps

    # --- Rearrange list from outer rings to inner rings --- #
    rings.out2in()
    hcomp.out2in()
    scomp.out2in()
    disk.out2in()

    # ----- Print ring information ----- #
    if display==1:
        # ----- Print (Start) ----- #
        print('--------------------------------------------------')
        print('R [Rg]: ', end='')
        for i_r, r in enumerate(rings.rs):
            if i_r==rings.n_ring-1:
                print('{:.1f}'.format(r))
            else:
                print('{:.1f}, '.format(r), end='')
        # ----- Print (End)   ----- #

        # ----- Print (Start) ----- #
        print('--------------------------------------------------')
        print('Disk ring [Rg]: ', end='')
        for i_r, r in enumerate(disk.rs):
            if i_r==len(disk.rs)-1:
                print('{:.1f}'.format(r))
            else:
                print('{:.1f}, '.format(r), end='')
        print('Soft Compton ring [Rg]: ', end='')
        for i_r, r in enumerate(scomp.rs):
            if i_r==len(scomp.rs)-1:
                print('{:.1f}'.format(r))
            else:
                print('{:.1f}, '.format(r), end='')
        print('Hard Compton ring [Rg]: ', end='')
        for i_r, r in enumerate(hcomp.rs):
            if i_r==len(hcomp.rs)-1:
                print('{:.1f}'.format(r))
            else:
                print('{:.1f}, '.format(r), end='')
        # ----- Print (End)   ----- #

        # ----- Print (Start) ----- #
        print('--------------------------------------------------')
        print('f_vis [Hz]: ', end='')
        for i_f, f_vis in enumerate(rings.fs_vis*bunit.c_rg):
            if i_f==rings.n_ring-1:
                print('{:.3f}'.format(f_vis))
            else:
                print('{:.3f}, '.format(f_vis), end='')
        # ----- Print (End)   ----- #

        # ----- Print (Start) ----- #
        print('--------------------------------------------------')
        print('v_rad [km/s]: ', end='')
        for i_v, v_rad in enumerate(rings.vs_rad*bunit.c):
            if i_v==rings.n_ring-1:
                print('{:.3f}'.format(v_rad))
            else:
                print('{:.3f}, '.format(v_rad), end='')
        # ----- Print (End)   ----- #

    # ------------------------------------------------------------------------------ #
    # ---------- PSD of mass accretion rate for each ring w/o propagation ---------- #
    # ------------------------------------------------------------------------------ #
    sigma=lf_var/np.sqrt(rings.n_dec) #\mu is fixed to unity. Probably this does not lose generality. (2021/07/26)
    flupro=FluPro()
    flupro.sigma=sigma
    flupro.f_set(fs=tifr.fs)
    flupro.f_vis_set(fs_vis=rings.fs_vis)
    flupro.psd_wo_prop(c_rg=bunit.c_rg)

    # ------------------------------------------------------------------------------- #
    # ---------- PSD of mass accretion rate for each ring with propagation ---------- #
    # ------------------------------------------------------------------------------- #
    flupro.dt=tifr.dt
    flupro.n_data=tifr.n_data
    flupro.psd_w_prop()
    if display==1:
        print('--------------------------------------------------')
        print('PSD of the mass accretion rate with propagation was successfully calculated.')

    # ------------------------------------------------------------------------------------------------------------- #
    # ---------- Calculation of w_flow=\int dE w(r_n, E), w_tot_flow := \sum _{r<r_ds} \int dE w(r_n, E) ---------- #
    # ------------------------------------------------------------------------------------------------------------- #
    ws_disk=np.zeros(len(disk.rs))
    ws_scomp=scomp.eps*spec.scomprr/scomp.eps_tot
    ws_hcomp=hcomp.eps*spec.hcomprr/hcomp.eps_tot
    ws=np.append(ws_disk, ws_scomp)
    ws=np.append(ws, ws_hcomp)
    # Integration over E
    # The energy range is not perfectly accurate.
    ws_flow=ws*(spec.e_maxrr-spec.e_minrr)  
    w_flow_tot=np.sum(ws_flow)

    ##############################################
    ########## Calculate power spectrum ##########
    ##############################################
    if quant==1:
        # ----------------------- #
        # ----- Energy band ----- #
        # ----------------------- #
        if display==1:
            print('--------------------------------------------------')
            print('Calculating power spectrum for')
            print('{0:.2f} - {1:.2f} keV'.format(spec.e_min, spec.e_max))
        md2fl=Mdot2Flux()
        md2fl.ene_set(e_min=spec.e_min, e_max=spec.e_max)
        #################################################################################
        ###### (2021/08/17) Preliminary (haphazard) prescription to set weight, ... #####
        ###### Smarter implementation will be performed.                            #####
        #################################################################################
        ws_disk =disk.eps *spec.disk /disk.eps_tot
        ws_scomp=scomp.eps*spec.scomp/scomp.eps_tot
        ws_hcomp=hcomp.eps*spec.hcomp/hcomp.eps_tot
        ws=np.append(ws_disk, ws_scomp)
        ws=np.append(ws, ws_hcomp)
        md2fl.ws=ws
        md2fl.w_tot=np.sum(md2fl.ws)

        md2fl.speceff_disk=spec.disk
        md2fl.speceff_sref=spec.sref
        md2fl.speceff_href=spec.href
        md2fl.mu_fl=spec.tot

        md2fl.psd_norm_set(dt=tifr.dt, n_data=tifr.n_data)
        f_dir=1. # fixed (not free parameter anymore and should be removed)
        f_rep=1.-f_dir # Fraction of the reprocessed component
        md2fl.norm_rep_set(f_rep=f_rep, dt0=dt0, w_flow_tot=w_flow_tot)

        lm2s_prop=flupro.psds_prop/flupro.norm_psd #|M_dot(r, f)|^2

        md2fl.psd_flux_rep_calc(fs=flupro.fs,\
                                n_r=rings.n_ring,\
                                ws_flow=ws_flow,\
                                lm2s=lm2s_prop,\
                                fs_vis=rings.fs_vis,\
                                dr_r=rings.dr_r,\
                                t0=t0,\
                                dt0=dt0,\
                                rg_c=bunit.rg_c)

        mus_fl=md2fl.mu_fl
        ws=md2fl.ws
        ws_tot=md2fl.w_tot
        norms_rep=md2fl.norm_rep
        psds_fl=md2fl.psd_fl

        if display==1:
            print('--------------------------------------------------')
            print('PSD of the flux was successfully calculated.')

        '''
        ###########################################
        ########## Simplest integratinon ##########
        ###########################################
        for i_f_data, f_data in enumerate(fs_data):
            i_f=np.abs(flupro.fs-f_data).argmin()
            psd=psds_fl[i_f]
            if i_f_data==0:
                psds=psd
            else:
                psds=np.append(psds, psd)

        n_f_data=len(fs_data)
        for i_f_data in range(n_f_data-1):
            # the most simple integration (area of trapezoid)
            flux=(psds[i_f_data]+psds[i_f_data+1])*(fs_data[i_f_data+1]-fs_data[i_f_data])/2.
            fluxes[i_f_data]=flux

        '''

        ################################################
        ########## More accurate integratinon ##########
        ################################################
        n_f_data=len(fs_data)
        for i_f_data in range(n_f_data-1):
            f_data_min=fs_data[i_f_data]
            f_data_max=fs_data[i_f_data+1]
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

        return flupro.fs, psds_fl

    ##############################################
    ########## Calculate cross spectrum ##########
    ##############################################
    elif quant in [2, 3, 4, 5, 6]:
        if display==1:
            print('--------------------------------------------------')
            print('Calculating cross spectrum for')
            print('{0:.2f} - {1:.2f} keV vs {2:.2f} - {3:.2f} keV'\
                  .format(spec.e_min, spec.e_max, spec.e_minr, spec.e_maxr))

        # ----- Weight ----- #
        md2fl=Mdot2Flux()

        ###################
        ### Energy band ###
        ###################
        md2fl.ene_set(e_min=spec.e_min, e_max=spec.e_max)
        #################################################################################
        ###### (2021/08/17) Preliminary (haphazard) prescription to set weight, ... #####
        ###### Smarter implementation will be performed.                            #####
        #################################################################################
        ws_disk =disk.eps *spec.disk /disk.eps_tot
        ws_scomp=scomp.eps*spec.scomp/scomp.eps_tot
        ws_hcomp=hcomp.eps*spec.hcomp/hcomp.eps_tot
        ws=np.append(ws_disk, ws_scomp)
        ws=np.append(ws, ws_hcomp)
        md2fl.ws=ws
        md2fl.w_tot=np.sum(md2fl.ws)
        md2fl.speceff_disk=spec.disk
        md2fl.speceff_sref=spec.sref
        md2fl.speceff_href=spec.href
        md2fl.mu_fl=spec.tot

        ######################
        ### Reference band ###
        ######################
        ws_disk =disk.eps *spec.diskr /disk.eps_tot
        ws_scomp=scomp.eps*spec.scompr/scomp.eps_tot
        ws_hcomp=hcomp.eps*spec.hcompr/hcomp.eps_tot
        ws=np.append(ws_disk, ws_scomp)
        ws=np.append(ws, ws_hcomp)
        md2fl.ws_ref=ws
        md2fl.w_tot_ref=np.sum(md2fl.ws_ref)
        md2fl.speceff_disk_ref=spec.diskr
        md2fl.speceff_sref_ref=spec.srefr
        md2fl.speceff_href_ref=spec.hrefr
        md2fl.mu_fl_ref=spec.totr

        # ----- Impulse response ----- #
        ### Channel-of-interest ###
        f_dir=1. # fixed (not free parameter anymore and should be removed)
        f_rep=1.-f_dir # Fraction of the reprocessed component
        md2fl.norm_rep_set(f_rep=f_rep, dt0=dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response

        ### Reference band ###
        md2fl.norm_rep_ref_set(f_rep=f_rep, dt0=dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response

        # ----- Normalization ----- #
        md2fl.csd_norm_set(dt=tifr.dt, n_data=tifr.n_data)

        # ----- CSD ----- #
        lm2s_prop=flupro.psds_prop/flupro.norm_psd #|M_dot(r, f)|^2

        md2fl.csd_flux_rep_calc(fs=flupro.fs,\
                                n_r=rings.n_ring,\
                                ws_rep=ws_flow,\
                                lm2s=lm2s_prop,\
                                fs_vis=rings.fs_vis,\
                                dr_r=rings.dr_r,\
                                t0=t0,\
                                dt0=dt0,\
                                rg_c=bunit.rg_c)

        # Re[CSD]
        if quant==2:
            csds_fl=np.real(md2fl.csd_fl)
        # Im[CSD]
        elif quant==3:
            csds_fl=np.imag(md2fl.csd_fl)
        # |CSD|
        elif quant==4:
            csds_fl=np.abs(md2fl.csd_fl)
        # Phase lag
        elif quant==5:
            md2fl.lagf_calc()
            csds_fl=md2fl.phi
        # Time lag
        elif quant==6:
            md2fl.lagf_calc()
            csds_fl=md2fl.tau

        if display==1:
            print('--------------------------------------------------')
            print('CSD of the flux was successfully calculated.')

        '''
        ###########################################
        ########## Simplest integratinon ##########
        ###########################################
        for i_f_data, f_data in enumerate(fs_data):
            i_f=np.abs(flupro.fs-f_data).argmin()
            csd=csds_fl[i_f]
            if i_f_data==0:
                csds=csd
            else:
                csds=np.append(csds, csd)

        n_f_data=len(fs_data)
        for i_f_data in range(n_f_data-1):
            # the most simple integration (area of trapezoid)
            flux=(csds[i_f_data]+csds[i_f_data+1])*(fs_data[i_f_data+1]-fs_data[i_f_data])/2.
            fluxes[i_f_data]=flux
        '''

        ################################################
        ########## More accurate integratinon ##########
        ################################################
        n_f_data=len(fs_data)
        for i_f_data in range(n_f_data-1):
            f_data_min=fs_data[i_f_data]
            f_data_max=fs_data[i_f_data+1]
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
