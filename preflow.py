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

    # Variability time-scale
    fs_vis_flow=f_vis_calc(r=rs_flow, lb=inpar.lb_flow, m=inpar.m_flow) #[c/Rg]
    fs_vis_disk=f_vis_calc(r=rs_disk, lb=inpar.lb_disk, m=inpar.m_disk) #[c/Rg]
    rings.fs_vis=np.append(fs_vis_flow, fs_vis_disk)
    # Accretion time-scale
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

    # --- Weight --- #
    # Energy band
    ws_hcomp=weight_calc_pivot(r=rs_hcomp, cc=inpar.cc_hcomp,\
        gamma=inpar.gamma_flow,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_mcomp=weight_calc_pivot(r=rs_mcomp, cc=inpar.cc_mcomp,\
        gamma=inpar.gamma_flow,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_scomp=weight_calc_pivot(r=rs_scomp, cc=inpar.cc_scomp,\
        gamma=inpar.gamma_flow,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_disk =weight_calc_pivot(r=rs_disk,  cc=inpar.cc_disk,\
        gamma=inpar.gamma_disk,   r_in=inpar.r_min,   stress=inpar.stress)
    ws=ws_disk
    ws=np.append(ws, ws_scomp)
    ws=np.append(ws, ws_mcomp)
    ws=np.append(ws, ws_hcomp)
    # Reference band
    ws_hcomp=weight_calc_pivot(r=rs_hcomp, cc=inpar.cc_hcompr,\
        gamma=inpar.gamma_flow,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_mcomp=weight_calc_pivot(r=rs_mcomp, cc=inpar.cc_mcompr,\
        gamma=inpar.gamma_flow,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_scomp=weight_calc_pivot(r=rs_scomp, cc=inpar.cc_scompr,\
        gamma=inpar.gamma_flow,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_disk =weight_calc_pivot(r=rs_disk,  cc=inpar.cc_diskr,\
        gamma=inpar.gamma_disk,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_r=ws_disk
    ws_r=np.append(ws_r, ws_scomp)
    ws_r=np.append(ws_r, ws_mcomp)
    ws_r=np.append(ws_r, ws_hcomp)
    # Reference band for reflection
    ws_hcomp=weight_calc_pivot(r=rs_hcomp, cc=inpar.cc_hcomprr,\
        gamma=inpar.gamma_flow,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_mcomp=weight_calc_pivot(r=rs_mcomp, cc=inpar.cc_mcomprr,\
        gamma=inpar.gamma_flow,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_scomp=weight_calc_pivot(r=rs_scomp, cc=inpar.cc_scomprr,\
        gamma=inpar.gamma_flow,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_disk =weight_calc_pivot(r=rs_disk,  cc=inpar.cc_diskrr,\
        gamma=inpar.gamma_disk,   r_in=inpar.r_min,   stress=inpar.stress)
    ws_rr=ws_disk
    ws_rr=np.append(ws_rr, ws_scomp)
    ws_rr=np.append(ws_rr, ws_mcomp)
    ws_rr=np.append(ws_rr, ws_hcomp)
    #print(ws)
    #print(ws_r)
    #print(ws_rr)

    # ----- Print ring information ----- #
    if inpar.display==1:
        print_ring_info(name='R [Rg]',                     xs=rings.rs,                digit=1)
        print_ring_info(name='Disk ring [Rg]',             xs=rs_disk,                 digit=1)
        print_ring_info(name='Soft Compton ring [Rg]',     xs=rs_scomp,                digit=1)
        print_ring_info(name='Mid Compton ring [Rg]',      xs=rs_mcomp,                digit=1)
        print_ring_info(name='Hard Compton ring [Rg]',     xs=rs_hcomp,                digit=1)
        print_ring_info(name='Variability frequency [Hz]', xs=rings.fs_vis*bunit.c_rg, digit=3)
        print_ring_info(name='Accretion frequency [Hz]',   xs=rings.fs_acc*bunit.c_rg, digit=3)
        print_ring_info(name='Radial velocity [km/s]',     xs=rings.vs_rad*bunit.c,    digit=3)

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
    #print(flupro.fs)
    #print(flupro.psds_prop[-1])

    # ------------------------------------------------------------------------------------------------------------- #
    # ---------- Calculation of w_flow=\int dE w(r_n, E), w_tot_flow := \sum _{r<r_ds} \int dE w(r_n, E) ---------- #
    # ------------------------------------------------------------------------------------------------------------- #
    #ws_disk=np.zeros(len(disk.rs))
    #ws_scomp=scomp.eps*spec.scomprr/scomp.eps_tot
    #ws_mcomp=mcomp.eps*spec.mcomprr/mcomp.eps_tot
    #ws_hcomp=hcomp.eps*spec.hcomprr/hcomp.eps_tot
    #ws=ws_disk
    #ws=np.append(ws, ws_scomp)
    #ws=np.append(ws, ws_mcomp)
    #ws=np.append(ws, ws_hcomp)
    ## Integration over E
    ## The energy range is not perfectly accurate.
    #ws_flow=ws*(spec.e_maxrr-spec.e_minrr)  
    #w_flow_tot=np.sum(ws_flow)

    # Integration over E
    #ws_flow=ws_rr*(inpar.e_maxrr-inpar.e_minrr)  
    #w_flow_tot=np.sum(ws_flow)

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
        #ws_disk =disk.eps *spec.disk /disk.eps_tot
        #ws_scomp=scomp.eps*spec.scomp/scomp.eps_tot
        #ws_mcomp=mcomp.eps*spec.mcomp/mcomp.eps_tot
        #ws_hcomp=hcomp.eps*spec.hcomp/hcomp.eps_tot
        #ws=ws_disk
        #ws=np.append(ws, ws_scomp)
        #ws=np.append(ws, ws_mcomp)
        #ws=np.append(ws, ws_hcomp)
        md2fl.ws=ws
        md2fl.w_tot=np.sum(md2fl.ws)

        #md2fl.speceff_disk=spec.disk
        #md2fl.speceff_sref=spec.sref
        #md2fl.speceff_mref=spec.mref
        #md2fl.speceff_href=spec.href
        #md2fl.mu_fl=spec.tot

        md2fl.speceff_disk=0.
        md2fl.speceff_sref=0.
        md2fl.speceff_mref=0.
        md2fl.speceff_href=0.
        md2fl.mu_fl=1. # <x(E,t)>=\sum _{r_n} w(r_n, E)=1

        md2fl.psd_norm_set(dt=tifr.dt, n_data=tifr.n_data)
        #f_dir=1. # fixed (not free parameter anymore and should be removed)
        #f_rep=1.-f_dir # Fraction of the reprocessed component
        #md2fl.norm_rep_set(f_rep=f_rep, dt0=dt0, w_flow_tot=w_flow_tot)
        #md2fl.norm_rep_set(dt0=inpar.dt0, w_flow_tot=w_flow_tot)

        md2fl.norm_rep=inpar.h0_rep # Normalization of the top-hat response function.

        lm2s_prop=flupro.psds_prop/flupro.norm_psd #|M_dot(r, f)|^2

        #md2fl.psd_flux_rep_calc(fs=flupro.fs,\
        #                        n_r=rings.n_ring,\
        #                        ws_flow=ws_rr,\
        #                        lm2s=lm2s_prop,\
        #                        fs_vis=rings.fs_vis,\
        #                        cds=rings.cds,\
        #                        xlag=inpar.xlag,\
        #                        dr_r=rings.dr_r,\
        #                        t0=inpar.t0,\
        #                        dt0=inpar.dt0,\
        #                        rg_c=bunit.rg_c)
        # fs_vis here is used to calculate accretion time!
        md2fl.psd_flux_rep_calc(fs=flupro.fs,\
                                n_r=rings.n_ring,\
                                ws_flow=ws_rr,\
                                lm2s=lm2s_prop,\
                                fs_vis=rings.fs_acc,\
                                cds=rings.cds,\
                                xlag=inpar.xlag,\
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

        #print(psds_fl)
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
        #ws_disk =disk.eps *spec.disk /disk.eps_tot
        #ws_scomp=scomp.eps*spec.scomp/scomp.eps_tot
        #ws_mcomp=mcomp.eps*spec.mcomp/mcomp.eps_tot
        #ws_hcomp=hcomp.eps*spec.hcomp/hcomp.eps_tot
        #ws=ws_disk
        #ws=np.append(ws, ws_scomp)
        #ws=np.append(ws, ws_mcomp)
        #ws=np.append(ws, ws_hcomp)
        #md2fl.ws=ws
        #md2fl.w_tot=np.sum(md2fl.ws)
        #md2fl.speceff_disk=spec.disk
        #md2fl.speceff_sref=spec.sref
        #md2fl.speceff_mref=spec.mref
        #md2fl.speceff_href=spec.href
        #md2fl.mu_fl=spec.tot

        md2fl.ws=ws
        md2fl.w_tot=np.sum(md2fl.ws)
        # Ignore reflection for now (2021/12/10)
        md2fl.speceff_disk=0.
        md2fl.speceff_sref=0.
        md2fl.speceff_mref=0.
        md2fl.speceff_href=0.
        md2fl.mu_fl=1. # <x(E,t)>=\sum _{r_n} w(r_n, E)=1

        ######################
        ### Reference band ###
        ######################
        #ws_disk =disk.eps *spec.diskr /disk.eps_tot
        #ws_scomp=scomp.eps*spec.scompr/scomp.eps_tot
        #ws_mcomp=mcomp.eps*spec.mcompr/mcomp.eps_tot
        #ws_hcomp=hcomp.eps*spec.hcompr/hcomp.eps_tot
        #ws=ws_disk
        #ws=np.append(ws, ws_scomp)
        #ws=np.append(ws, ws_mcomp)
        #ws=np.append(ws, ws_hcomp)
        #md2fl.ws_ref=ws
        #md2fl.w_tot_ref=np.sum(md2fl.ws_ref)
        #md2fl.speceff_disk_ref=spec.diskr
        #md2fl.speceff_sref_ref=spec.srefr
        #md2fl.speceff_mref_ref=spec.mrefr
        #md2fl.speceff_href_ref=spec.hrefr
        #md2fl.mu_fl_ref=spec.totr

        md2fl.ws_ref=ws_r
        md2fl.w_tot_ref=np.sum(md2fl.ws_ref)
        # Ignore reflection for now (2021/12/10)
        md2fl.speceff_disk_ref=0.
        md2fl.speceff_sref_ref=0.
        md2fl.speceff_mref_ref=0.
        md2fl.speceff_href_ref=0.
        md2fl.mu_fl_ref=1. # <x(E,t)>=\sum _{r_n} w(r_n, E)=1

        # ----- Impulse response ----- #
        ### Channel-of-interest ###
        #f_dir=1. # fixed (not free parameter anymore and should be removed)
        #f_rep=1.-f_dir # Fraction of the reprocessed component
        #md2fl.norm_rep_set(f_rep=f_rep, dt0=dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response
        #md2fl.norm_rep_set(dt0=inpar.dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response
        md2fl.norm_rep=inpar.h0_rep # Normalization of the top-hat response function.

        ### Reference band ###
        #md2fl.norm_rep_ref_set(f_rep=f_rep, dt0=dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response
        #md2fl.norm_rep_ref_set(dt0=inpar.dt0, w_flow_tot=w_flow_tot) # C(E) in the impulse response
        md2fl.norm_rep_ref=inpar.h0_repr # Normalization of the top-hat response function.

        # ----- Normalization ----- #
        md2fl.csd_norm_set(dt=tifr.dt, n_data=tifr.n_data)

        # ----- CSD ----- #
        lm2s_prop=flupro.psds_prop/flupro.norm_psd #|M_dot(r, f)|^2

        #md2fl.csd_flux_rep_calc(fs=flupro.fs,\
        #                        n_r=rings.n_ring,\
        #                        ws_rep=ws_rr,\
        #                        lm2s=lm2s_prop,\
        #                        fs_vis=rings.fs_vis,\
        #                        cds=rings.cds,\
        #                        xlag=inpar.xlag,\
        #                        dr_r=rings.dr_r,\
        #                        t0=inpar.t0,\
        #                        dt0=inpar.dt0,\
        #                        rg_c=bunit.rg_c)
        # fs_vis here is used to calculate accretion time!
        md2fl.csd_flux_rep_calc(fs=flupro.fs,\
                                n_r=rings.n_ring,\
                                ws_rep=ws_rr,\
                                lm2s=lm2s_prop,\
                                fs_vis=rings.fs_acc,\
                                cds=rings.cds,\
                                xlag=inpar.xlag,\
                                dr_r=rings.dr_r,\
                                t0=inpar.t0,\
                                dt0=inpar.dt0,\
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
