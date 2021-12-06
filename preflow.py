from preflow_h import *
from scipy import integrate

def preflow(engs, params, fluxes):
    # ------------------------------------- #
    # ---------- Initial setting ---------- #
    # ------------------------------------- #
    inpar=SetParameter()
    inpar.set_inpar(pars=params, es=engs)
    inpar.check_validity()
    
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

    # ------------------------------- #
    # ---------- Set flux  ---------- #
    # ------------------------------- #
    spec=FluxData()
    spec.set_flux(e_min  =inpar.e_min,\
                  e_max  =inpar.e_max,\
                  disk   =inpar.frac_disk,\
                  scomp  =inpar.frac_scomp,\
                  mcomp  =inpar.frac_mcomp,\
                  hcomp  =inpar.frac_hcomp,\
                  sref   =inpar.frac_sref,\
                  mref   =inpar.frac_mref,\
                  href   =inpar.frac_href,\
                  tot    =inpar.frac_tot,\
                  e_minr =inpar.e_minr,\
                  e_maxr =inpar.e_maxr,\
                  diskr  =inpar.frac_diskr,\
                  scompr =inpar.frac_scompr,\
                  mcompr =inpar.frac_mcompr,\
                  hcompr =inpar.frac_hcompr,\
                  srefr  =inpar.frac_srefr,\
                  mrefr  =inpar.frac_mrefr,\
                  hrefr  =inpar.frac_hrefr,\
                  totr   =inpar.frac_totr,\
                  e_minrr=inpar.e_minrr,\
                  e_maxrr=inpar.e_maxrr,\
                  scomprr=inpar.frac_scomprr,\
                  mcomprr=inpar.frac_mcomprr,\
                  hcomprr=inpar.frac_hcomprr)

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

    # ------------------------------------------ #
    # ---------- Assign flux to ring  ---------- #
    # ------------------------------------------ #
    r_min_eps=inpar.r_in
    ###########################
    ### Hard Comptonization ###
    ###########################
    hcomp=Ring2Spec()
    hcomp.r_min=inpar.r_in
    hcomp.r_max=inpar.r_mh
    hcomp.lb=inpar.lb_flow
    hcomp.m=inpar.m_flow
    hcomp.epsilon_par_set(stress=inpar.stress, gamma=inpar.gamma, r_min=r_min_eps)
    hcomp.ring_assign(rs=rings.rs,\
                      rs_min=rings.rs_min,\
                      wids=rings.wids,\
                      drs=rings.drs)
    hcomp.f_vis_set()
    hcomp.v_rad_set()
    hcomp.epsilon_set()

    ##########################
    ### Mid Comptonization ###
    ##########################
    mcomp=Ring2Spec()
    mcomp.r_min=inpar.r_mh
    mcomp.r_max=inpar.r_sm
    mcomp.lb=inpar.lb_flow
    mcomp.m=inpar.m_flow
    mcomp.epsilon_par_set(stress=inpar.stress, gamma=inpar.gamma, r_min=r_min_eps)
    mcomp.ring_assign(rs=rings.rs,\
                      rs_min=rings.rs_min,\
                      wids=rings.wids,\
                      drs=rings.drs)
    mcomp.f_vis_set()
    mcomp.v_rad_set()
    mcomp.epsilon_set()

    ###########################
    ### Soft Comptonization ###
    ###########################
    scomp=Ring2Spec()
    scomp.r_min=inpar.r_sm
    scomp.r_max=inpar.r_ds
    scomp.lb=inpar.lb_flow
    scomp.m=inpar.m_flow
    scomp.epsilon_par_set(stress=inpar.stress, gamma=inpar.gamma, r_min=r_min_eps)
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
    disk.r_min=inpar.r_ds
    disk.r_max=inpar.r_out
    disk.lb=inpar.lb_disk
    disk.m=inpar.m_disk
    disk.epsilon_par_set(stress=inpar.stress, gamma=inpar.gamma, r_min=r_min_eps)
    disk.ring_assign(rs=rings.rs,\
                     rs_min=rings.rs_min,\
                     wids=rings.wids,\
                     drs=rings.drs)

    # -------------------------------------- #
    # ---------- Original (begin) ---------- #
    # -------------------------------------- #
    disk.f_vis_set()
    # -------------------------------------- #
    # ---------- Original (end)   ---------- #
    # -------------------------------------- #
    disk.v_rad_set()
    disk.epsilon_set()

    # --- Give info on viscous frequency and radial velocity to Flow2Ring class instance --- #
    fs_vis=hcomp.fs_vis
    fs_vis=np.append(fs_vis, mcomp.fs_vis)
    fs_vis=np.append(fs_vis, scomp.fs_vis)
    fs_vis=np.append(fs_vis, disk.fs_vis)
    rings.fs_vis=fs_vis

    vs_rad=hcomp.vs_rad
    vs_rad=np.append(vs_rad, mcomp.vs_rad)
    vs_rad=np.append(vs_rad, scomp.vs_rad)
    vs_rad=np.append(vs_rad, disk.vs_rad)
    rings.vs_rad=vs_rad

    eps=hcomp.eps
    eps=np.append(eps, mcomp.eps)
    eps=np.append(eps, scomp.eps)
    eps=np.append(eps, disk.eps)
    rings.eps=eps

    # --- Rearrange list from outer rings to inner rings --- #
    rings.out2in()
    hcomp.out2in()
    scomp.out2in()
    disk.out2in()

    # ----- Print ring information ----- #
    if inpar.display==1:
        print_ring_info(name='R [Rg]',                 xs=rings.rs,                digit=1)
        print_ring_info(name='Disk ring [Rg]',         xs=disk.rs,                 digit=1)
        print_ring_info(name='Soft Compton ring [Rg]', xs=scomp.rs,                digit=1)
        print_ring_info(name='Mid Compton ring [Rg]',  xs=mcomp.rs,                digit=1)
        print_ring_info(name='Hard Compton ring [Rg]', xs=hcomp.rs,                digit=1)
        print_ring_info(name='Viscous frequency [Hz]', xs=rings.fs_vis*bunit.c_rg, digit=3)
        print_ring_info(name='Radial velocity [km/s]', xs=rings.vs_rad*bunit.c,    digit=3)
        
    # ------------------------------------------------------------------------------ #
    # ---------- PSD of mass accretion rate for each ring w/o propagation ---------- #
    # ------------------------------------------------------------------------------ #
    sigma=inpar.lf_var/np.sqrt(rings.n_dec) #\mu is fixed to unity. Probably this does not lose generality. (2021/07/26)
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
    if inpar.display==1:
        print('--------------------------------------------------')
        print('PSD of the mass accretion rate with propagation was successfully calculated.')

    # ------------------------------------------------------------------------------------------------------------- #
    # ---------- Calculation of w_flow=\int dE w(r_n, E), w_tot_flow := \sum _{r<r_ds} \int dE w(r_n, E) ---------- #
    # ------------------------------------------------------------------------------------------------------------- #
    ws_disk=np.zeros(len(disk.rs))
    ws_scomp=scomp.eps*spec.scomprr/scomp.eps_tot
    ws_mcomp=mcomp.eps*spec.mcomprr/mcomp.eps_tot
    ws_hcomp=hcomp.eps*spec.hcomprr/hcomp.eps_tot
    ws=ws_disk
    ws=np.append(ws, ws_scomp)
    ws=np.append(ws, ws_mcomp)
    ws=np.append(ws, ws_hcomp)
    # Integration over E
    # The energy range is not perfectly accurate.
    ws_flow=ws*(spec.e_maxrr-spec.e_minrr)  
    w_flow_tot=np.sum(ws_flow)

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
            print('{0:.2f} - {1:.2f} keV'.format(spec.e_min, spec.e_max))
        md2fl=Mdot2Flux()
        md2fl.ene_set(e_min=spec.e_min, e_max=spec.e_max)
        #################################################################################
        ###### (2021/08/17) Preliminary (haphazard) prescription to set weight, ... #####
        ###### Smarter implementation will be performed.                            #####
        #################################################################################
        ws_disk =disk.eps *spec.disk /disk.eps_tot
        ws_scomp=scomp.eps*spec.scomp/scomp.eps_tot
        ws_mcomp=mcomp.eps*spec.mcomp/mcomp.eps_tot
        ws_hcomp=hcomp.eps*spec.hcomp/hcomp.eps_tot
        ws=ws_disk
        ws=np.append(ws, ws_scomp)
        ws=np.append(ws, ws_mcomp)
        ws=np.append(ws, ws_hcomp)
        md2fl.ws=ws
        md2fl.w_tot=np.sum(md2fl.ws)

        md2fl.speceff_disk=spec.disk
        md2fl.speceff_sref=spec.sref
        md2fl.speceff_mref=spec.mref
        md2fl.speceff_href=spec.href
        md2fl.mu_fl=spec.tot

        md2fl.psd_norm_set(dt=tifr.dt, n_data=tifr.n_data)
        #f_dir=1. # fixed (not free parameter anymore and should be removed)
        #f_rep=1.-f_dir # Fraction of the reprocessed component
        #md2fl.norm_rep_set(f_rep=f_rep, dt0=dt0, w_flow_tot=w_flow_tot)
        md2fl.norm_rep_set(dt0=inpar.dt0, w_flow_tot=w_flow_tot)

        lm2s_prop=flupro.psds_prop/flupro.norm_psd #|M_dot(r, f)|^2

        md2fl.psd_flux_rep_calc(fs=flupro.fs,\
                                n_r=rings.n_ring,\
                                ws_flow=ws_flow,\
                                lm2s=lm2s_prop,\
                                fs_vis=rings.fs_vis,\
                                dr_r=rings.dr_r,\
                                t0=inpar.t0,\
                                dt0=inpar.dt0,\
                                rg_c=bunit.rg_c)

        #mus_fl=md2fl.mu_fl
        #ws=md2fl.ws
        #ws_tot=md2fl.w_tot
        #norms_rep=md2fl.norm_rep
        psds_fl=md2fl.psd_fl

        if inpar.display==1:
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

        return flupro.fs, psds_fl

    ##############################################
    ########## Calculate cross spectrum ##########
    ##############################################
    elif inpar.quant in [2, 3, 4, 5, 6]:
        if inpar.display==1:
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
        ws_mcomp=mcomp.eps*spec.mcomp/mcomp.eps_tot
        ws_hcomp=hcomp.eps*spec.hcomp/hcomp.eps_tot
        ws=ws_disk
        ws=np.append(ws, ws_scomp)
        ws=np.append(ws, ws_mcomp)
        ws=np.append(ws, ws_hcomp)
        md2fl.ws=ws
        md2fl.w_tot=np.sum(md2fl.ws)
        md2fl.speceff_disk=spec.disk
        md2fl.speceff_sref=spec.sref
        md2fl.speceff_mref=spec.mref
        md2fl.speceff_href=spec.href
        md2fl.mu_fl=spec.tot

        ######################
        ### Reference band ###
        ######################
        ws_disk =disk.eps *spec.diskr /disk.eps_tot
        ws_scomp=scomp.eps*spec.scompr/scomp.eps_tot
        ws_mcomp=mcomp.eps*spec.mcompr/mcomp.eps_tot
        ws_hcomp=hcomp.eps*spec.hcompr/hcomp.eps_tot
        ws=ws_disk
        ws=np.append(ws, ws_scomp)
        ws=np.append(ws, ws_mcomp)
        ws=np.append(ws, ws_hcomp)
        md2fl.ws_ref=ws
        md2fl.w_tot_ref=np.sum(md2fl.ws_ref)
        md2fl.speceff_disk_ref=spec.diskr
        md2fl.speceff_sref_ref=spec.srefr
        md2fl.speceff_mref_ref=spec.mrefr
        md2fl.speceff_href_ref=spec.hrefr
        md2fl.mu_fl_ref=spec.totr

        # ----- Impulse response ----- #
        ### Channel-of-interest ###
        #f_dir=1. # fixed (not free parameter anymore and should be removed)
        #f_rep=1.-f_dir # Fraction of the reprocessed component
        #md2fl.norm_rep_set(f_rep=f_rep, dt0=dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response
        md2fl.norm_rep_set(dt0=inpar.dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response

        ### Reference band ###
        #md2fl.norm_rep_ref_set(f_rep=f_rep, dt0=dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response
        md2fl.norm_rep_ref_set(dt0=inpar.dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response

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
                                t0=inpar.t0,\
                                dt0=inpar.dt0,\
                                rg_c=bunit.rg_c)

        # Re[CSD]
        if inpar.quant==2:
            csds_fl=np.real(md2fl.csd_fl)
        # Im[CSD]
        elif inpar.quant==3:
            csds_fl=np.imag(md2fl.csd_fl)
        # |CSD|
        elif inpar.quant==4:
            csds_fl=np.abs(md2fl.csd_fl)
        # Phase lag
        elif inpar.quant==5:
            md2fl.lagf_calc()
            csds_fl=md2fl.phi
        # Time lag
        elif inpar.quant==6:
            md2fl.lagf_calc()
            csds_fl=md2fl.tau

        if inpar.display==1:
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
