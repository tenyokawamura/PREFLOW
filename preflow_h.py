import sys
import cmath
from fourier_h import *

# -------------------------- #
# ----- Read parameter ----- #
# -------------------------- #
class SetParameter:
    def __init__(self):
        self.set_inpar_done=False
        self.check_validity_done=False

    def set_inpar(self, pars, es):
        self.mass        =pars[0]  # BH mass [solar mass]
        self.r_in        =pars[1]  # Inner radius of hard Compton (inner radius of hot flow) [Rg]
        self.r_mh        =pars[2]  # Radius at which viscous frequency prescriptions changes [Rg]
        self.r_sm        =pars[3]  # Transition radius between hard and soft Compton [Rg]
        self.r_ds        =pars[4]  # Transition radius between soft Compton and variable disk [Rg]
        self.r_out       =pars[5]  # Transition radius between soft Compton and variable disk (outer radius of hot flow) [Rg]
        self.n_ring      =pars[6]  # Outer radius of variable disk [Rg]
        self.tref        =pars[7]  # Start time of reflection impulse response [sec]
        self.dtref       =pars[8]  # Time width of reflection impulse response [sec]
        self.lf_var      =pars[9]  # Fractional variability of mass accretion rate in radial decade [-]
        self.lb_disk     =pars[10] # B_{disk} [-]
        self.m_disk      =pars[11] # m_{disk} [-]
        self.lb_flow     =pars[12] # B_{flow} [-]
        self.m_flow      =pars[13] # m_{flow} [-]
        self.stress      =pars[14] # 1: stressed, 2: stress-free in emissivity
        self.gamma       =pars[15] # Radial index of emissivity [-]
        self.e_min       =pars[16] # Lower bound of energy band [keV] (unused)
        self.e_max       =pars[17] # Upper bound of energy band [keV] (unused)
        self.frac_disk   =pars[18] # Fraction of variable disk in the energy band [counts keV^-1 s^-1]
        self.frac_scomp  =pars[19] # Fraction of soft Compton in the energy band [counts keV^-1 s^-1]
        self.frac_mcomp  =pars[20] # Fraction of mid Compton in the energy band [counts keV^-1 s^-1]
        self.frac_hcomp  =pars[21] # Fraction of hard Compton in the energy band [counts keV^-1 s^-1]
        self.frac_sref   =pars[22] # Fraction of soft reflection in the energy band [counts keV^-1 s^-1]
        self.frac_mref   =pars[23] # Fraction of mid reflection in the energy band [counts keV^-1 s^-1]
        self.frac_href   =pars[24] # Fraction of hard reflection in the energy band [counts keV^-1 s^-1]
        self.e_minr      =pars[25] # Lower bound of reference band [keV] (unused)
        self.e_maxr      =pars[26] # Upper bound of reference band [keV] (unused)
        self.frac_diskr  =pars[27] # Fraction of variable disk in the reference band [counts keV^-1 s^-1]
        self.frac_scompr =pars[28] # Fraction of soft Compton in the reference band [counts keV^-1 s^-1]
        self.frac_mcompr =pars[29] # Fraction of mid Compton in the energy band [counts keV^-1 s^-1]
        self.frac_hcompr =pars[30] # Fraction of hard Compton in the reference band [counts keV^-1 s^-1]
        self.frac_srefr  =pars[31] # Fraction of soft reflection in the reference band [counts keV^-1 s^-1]
        self.frac_mrefr  =pars[32] # Fraction of mid reflection in the reference band [counts keV^-1 s^-1]
        self.frac_hrefr  =pars[33] # Fraction of hard reflection in the reference band [counts keV^-1 s^-1]
        self.e_minrr     =pars[34] # Lower bound of reference band 'for reflection' [keV] (unused)
        self.e_maxrr     =pars[35] # Upper bound of reference band 'for reflection' [keV] (unused)
        self.frac_scomprr=pars[36] # Soft Compton in the reference band 'for reflection' [counts kev^-1 s^-1]
        self.frac_mcomprr=pars[37] # Mid Compton in the reference band 'for reflection' [counts kev^-1 s^-1]
        self.quant       =pars[38]
            # 1: power spectrum 
            # 2: real part of cross spectrum
            # 3: imaginary part of cross spectrum
            # 4: absolute value of cross spectrum
            # 5: phase lag (Positive lag means reference band lagging behind energy band.)
            # 6: time lag  (Positive lag means reference band lagging behind energy band.)
        self.display     =pars[39] # 1: display, 2: not display

        # PREFLOW model is a timing model!
        # Energy in XSPEC corresponds to Fourier frequency in preflow.
        self.fs_data=es

        ### Parameters, which are no longer free. ###
        # Total fraction in the energy band [counts keV^-1 s^-1], i.e., 1.
        self.frac_tot =1. 
        # Total fraction in the reference band [counts keV^-1 s^-1], 1.e., 1.
        self.frac_totr=1. 
        # Fraction of hard Compton in the hot flow in the reference band 
        #'for reflection' [counts kev^-1 s^-1]
        self.frac_hcomprr=1.-(self.frac_scomprr+self.frac_mcomprr) 

        # Impulse response of reflection
        self.t0=self.tref+(self.dtref/2.)
        self.dt0=self.dtref

        self.set_inpar_done=True

    def check_validity(self):
        '''
        #######################################
        ### Need to be updated (2021/12/06) ###
        #######################################
        ### Flux ###
        # Energy band
        if frac_disk+frac_scomp+frac_mcomp+frac_hcomp+frac_sref+frac_href<=frac_tot:
            pass
        elif frac_disk>frac_tot:
            print('Warning: frac_disk>1 --> frac_disk=1, frac_scomp=frac_hcomp=frac_sref=frac_href=0')
            frac_disk=frac_tot
            frac_scomp=0
            frac_hcomp=0
            frac_sref=0
            frac_href=0
        elif frac_disk+frac_scomp>frac_tot:
            print('Warning: frac_disk+frac_scomp>1 --> frac_scomp=1-frac_disk, frac_hcomp=frac_sref=frac_href=0')
            frac_scomp=frac_tot-frac_disk
            frac_hcomp=0
            frac_sref=0
            frac_href=0
        elif frac_disk+frac_scomp+frac_hcomp>frac_tot:
            print('Warning: frac_disk+frac_scomp+frac_hcomp>1 --> frac_hcomp=1-(frac_disk+frac_scomp), frac_sref=frac_href=0')
            frac_hcomp=frac_tot-(frac_disk+frac_scomp)
            frac_sref=0
            frac_href=0
        elif frac_disk+frac_scomp+frac_hcomp+frac_sref>frac_tot:
            print('Warning: frac_disk+frac_scomp+frac_hcomp+frac_sref>1 --> frac_sref=1-(frac_disk+frac_scomp+frac_hcomp), frac_href=0')
            frac_sref=frac_tot-(frac_disk+frac_scomp+frac_hcomp)
            frac_href=0
        elif frac_disk+frac_scomp+frac_hcomp+frac_sref+frac_href>frac_tot:
            print('Warning: frac_disk+frac_scomp+frac_hcomp+frac_sref>1 --> frac_href=1-(frac_disk+frac_scomp+frac_hcomp+frac_sref)')
            frac_href=frac_tot-(frac_disk+frac_scomp+frac_hcomp+frac_sref)
        else:
            print('Error')

        # Reference band
        if frac_diskr+frac_scompr+frac_hcompr+frac_srefr+frac_hrefr<=frac_tot:
            pass
        elif frac_diskr>frac_tot:
            print('Warning: frac_diskr>1 --> frac_diskr=1, frac_scompr=frac_hcompr=frac_srefr=frac_hrefr=0')
            frac_diskr=frac_tot
            frac_scompr=0.
            frac_hcompr=0.
            frac_srefr=0.
            frac_hrefr=0.
        elif frac_diskr+frac_scompr>frac_tot:
            print('Warning: frac_diskr+frac_scompr>1 --> frac_scompr=1-frac_diskr, frac_hcompr=frac_srefr=frac_hrefr=0')
            frac_scompr=frac_tot-frac_diskr
            frac_hcompr=0.
            frac_srefr=0.
            frac_hrefr=0.
        elif frac_diskr+frac_scompr+frac_hcompr>frac_tot:
            print('Warning: frac_diskr+frac_scompr+frac_hcompr>1 --> frac_hcompr=1-(frac_diskr+frac_scompr), frac_srefr=frac_hrefr=0')
            frac_hcompr=frac_tot-(frac_diskr+frac_scompr)
            frac_srefr=0.
            frac_hrefr=0.
        elif frac_diskr+frac_scompr+frac_hcompr+frac_srefr>frac_tot:
            print('Warning: frac_diskr+frac_scompr+frac_hcompr+frac_srefr>1 --> frac_srefr=1-(frac_diskr+frac_scompr+frac_hcompr), frac_hrefr=0')
            frac_srefr=frac_tot-(frac_diskr+frac_scompr+frac_hcompr)
            frac_hrefr=0.
        elif frac_diskr+frac_scompr+frac_hcompr+frac_srefr+frac_hrefr>frac_tot:
            print('Warning: frac_diskr+frac_scompr+frac_hcompr+frac_srefr>1 --> frac_hrefr=1-(frac_diskr+frac_scompr+frac_hcompr+frac_srefr)')
            frac_hrefr=frac_tot-(frac_diskr+frac_scompr+frac_hcompr+frac_srefr)
        else:
            print('Error')

        # Reference band 'for reflection'
        if frac_scomprr<=1.:
            pass
        else:
            print('Warning: frac_scomprr>1 --> frac_scomprr=1, frac_hcomprr=0')
            frac_scomprr=1.
            frac_hcomprr=0.
        '''

        ### Geometry ###
        self.r_in, self.r_mh =check_validity_geo(\
            r_in=self.r_in, r_out=self.r_mh,\
            name_r_in='rin', name_r_out='rmh')
        self.r_mh, self.r_sm =check_validity_geo(\
            r_in=self.r_mh, r_out=self.r_sm,\
            name_r_in='rmh', name_r_out='rsm')
        self.r_sm, self.r_ds =check_validity_geo(\
            r_in=self.r_sm, r_out=self.r_ds,\
            name_r_in='rsm', name_r_out='rds')
        self.r_ds, self.r_out=check_validity_geo(\
            r_in=self.r_ds, r_out=self.r_out,\
            name_r_in='rds', name_r_out='rout')

        self.check_validity_done=True

def check_validity_geo(r_in, r_out, name_r_in, name_r_out):
    if r_in>r_out:
        print('Warning: {0}>{1} --> {0}={1}'.format(name_r_in, name_r_out))
        r_out=r_in
    return r_in, r_out

# -------------------------- #
# ----- Set basic unit ----- #
# -------------------------- #
class BasicUnit:
    def __init__(self):
        self.set_unit_done=False

    def set_unit(self, mass):
        # mass is in units of solar mass
        self.m=mass              # M [solar mass]
        self.c=3.e5              # c [km/s]
        self.rg=1.5*self.m       # Rg [km]
        self.rg_c=self.rg/self.c # Rg/c [s]
        self.c_rg=self.c/self.rg # c/Rg [Hz]
        self.set_unit_done=True

# ------------------------------ #
# ----- Set time/frequency ----- #
# ------------------------------ #
class TimeFrequency:
    def __init__(self):
        self.n_data=0
        self.dt=0

    def set_par(self, f_data_min, f_data_max):
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
        self.dt=dt

        # ----- Decide n_data ----- #
        n_data=1
        while True:
            n_data*=2
            f_min=1./(n_data*dt)
            if f_min<f_data_min:
                break
        self.n_data=n_data

    def t_set(self):
        self.lt=self.n_data*self.dt #[s]
        self.ts=np.arange(0, self.lt, self.dt)
        if len(self.ts)!=self.n_data:
            print('Error')
            sys.exit()

    def f_set(self):
        self.f_min=1./self.lt      #[Hz]
        self.f_max=1./(2.*self.dt) #[Hz]
        self.df=self.f_min #[Hz]
        self.fs=np.arange(self.f_min, self.f_max+self.df, self.df)
        if (2*len(self.fs)+1 - self.n_data) in [0, 1] == False:
            print('Error')
            print(2*len(self.fs)+1, self.n_data)
            sys.exit()

# ------------------------- #
# ----- Spectral data ----- #
# ------------------------- #
class FluxData:
    def __init__(self):
        self.set_flux_done=False

    def set_flux(self,\
                 e_min,   e_max,   disk,    scomp,   mcomp,   hcomp,   sref,  mref,  href,  tot,\
                 e_minr,  e_maxr,  diskr,   scompr,  mcompr,  hcompr,  srefr, mrefr, hrefr, totr,\
                 e_minrr, e_maxrr, scomprr, mcomprr, hcomprr):
        #################################################################################
        ##### Unit: [counts keV^-1 s^-1]                                            #####
        ##### But this unit is not intuitive and it is hard to give initial values. #####
        ##### This will be modified for easy use.                                   #####
        #################################################################################
        # ----- Energy band ----- #
        self.e_min=e_min
        self.e_max=e_max
        self.disk=disk
        self.scomp=scomp
        self.mcomp=mcomp
        self.hcomp=hcomp
        self.sref=sref
        self.mref=mref
        self.href=href
        self.tot=tot

        # ----- Reference band ----- #
        self.e_minr=e_minr
        self.e_maxr=e_maxr
        self.diskr=diskr
        self.scompr=scompr
        self.mcompr=mcompr
        self.hcompr=hcompr
        self.srefr=srefr
        self.mrefr=mrefr
        self.hrefr=hrefr
        self.totr=totr

        # ----- Reference band 'for refection' ----- #
        self.e_minrr=e_minrr
        self.e_maxrr=e_maxrr
        self.scomprr=scomprr
        self.mcomprr=mcomprr
        self.hcomprr=hcomprr

        self.set_flux_done=True

# ------------------------------------ #
# ----- Separate flow into rings ----- #
# ------------------------------------ #
class Flow2Ring:
    def __init__(self):
        self.n_ring=0
        self.r_in=0
        self.r_out=0
        self.k=0
        self.n_dec=0
        self.set=False

    # Initial setting
    def par_set(self):
        self.rs_min=np.zeros(self.n_ring)
        self.rs=np.zeros(self.n_ring)
        self.wids=np.zeros(self.n_ring)
        self.drs=np.zeros(self.n_ring-1)
        self.fs_vis=np.zeros(self.n_ring)
        self.vs_rad=np.zeros(self.n_ring)
        self.eps=np.zeros(self.n_ring)
        self.set=True

    # Coefficient for logarithmical separation
    def coeff_calc(self):
        if self.r_out==0:
            print('Error: geometry is not set.')
            sys.exit()
        k_min=0
        k_max=10.
        dk=2.e-4
        ks=np.arange(k_min, k_max, dk)
        for i, k in enumerate(ks):
            r_min_out=self.r_in*(1.+k)**(self.n_ring-1)
            if r_min_out>self.r_out:
                self.k=ks[i-1]
                break
            elif i==len(ks)-1:
                print('Error: no appropriate coefficient is found.')
                sys.exit()

    # Lower boundary of ring
    def lowboundr_calc(self):
        if self.set==False:
            print('Error: parameter is not set.')
            sys.exit()
        if self.k==0:
            print('Error: k=0')
            sys.exit()
        for i in range(self.n_ring):
            if i==0:
                self.rs_min[i]=self.r_in
            else:
                self.rs_min[i]=self.rs_min[i-1]*(1.+self.k)

    # Center position of ring
    def centerr_calc(self):
        if self.set==False:
            print('Error: parameter is not set.')
            sys.exit()
        for i in range(self.n_ring):
            if i==self.n_ring-1:
                self.rs[i]=(self.rs_min[i]+self.r_out)/2.
            else:
                self.rs[i]=(self.rs_min[i]+self.rs_min[i+1])/2.

    # Width of ring
    def width_calc(self):
        if self.set==False:
            print('Error: parameter is not set.')
            sys.exit()
        for i in range(self.n_ring):
            if i==self.n_ring-1:
                self.wids[i]=self.r_out-self.rs_min[i]
            else:
                self.wids[i]=self.rs_min[i+1]-self.rs_min[i]

    # Interval of ring (between centers)
    def interval_calc(self):
        if self.set==False:
            print('Error: parameter is not set.')
            sys.exit()
        self.drs=np.roll(self.rs, -1)-self.rs
        self.drs=self.drs[:-1]
        self.dr_r=self.drs[0]/self.rs[0]

    # Calculate Ndec (Number of rings per decade)
    def n_ring_dec_calc(self):
        if self.r_in==0:
            print('Error: geometry is not set.')
            sys.exit()
        i=0
        while True:
            r_min=self.r_in*(1.+self.k)**(i)
            if r_min>10.*self.r_in:
                self.n_dec=i+1
                break
            else:
                i+=1

    # Final operation!
    def out2in(self):
        self.rs=self.rs[::-1]
        self.rs_min=self.rs_min[::-1]
        self.wids=self.wids[::-1]
        self.drs=self.drs[::-1]
        self.fs_vis=self.fs_vis[::-1]
        self.vs_rad=self.vs_rad[::-1]
        self.eps=self.eps[::-1]

# ----------------------------------- #
# ----- Assign spectrum to ring ----- #
# ----------------------------------- #
class Ring2Spec:
    def __init__(self):
        self.r_min=0
        self.r_max=0
        self.lb=0
        self.m=0
        self.ring_assign_done=False
        self.f_vis_set_done=False
        self.epsilon_par_set_done=False

    # Parameter for emissivity
    def epsilon_par_set(self, stress, gamma, r_min):
        if stress==1:
            self.stress=True  # stressed
        elif stress==2:
            self.stress=False # stress-free
        self.gamma_eps=gamma
        self.r_min_eps=r_min
        self.epsilon_par_set_done=True

    # Assign ring to spectrum
    def ring_assign(self, rs, rs_min, wids, drs):
        if self.r_min==0 and self.r_max==0:
            print('Error: geometry is not set.')
            sys.exit()
        start=False
        end=False
        n_r=len(rs)
        for i in range(n_r):
            if self.r_min<=rs_min[i] and start==False:
                i_start=i
                start=True

            if self.r_max<rs_min[i] and end==False:
                i_end=i
                end=True
                break

            if i==n_r-1:
                i_end=n_r

        self.rs=rs[i_start:i_end]
        self.rs_min=rs_min[i_start:i_end]
        self.wids=wids[i_start:i_end]
        self.drs=drs[i_start:i_end]
        self.ring_assign_done=True

    def f_vis_set(self):
        #if self.lb==0 or self.c_rg==0:
        #    print('Error: parameter is not set.')
        #    sys.exit()
        if self.ring_assign_done==False:
            print('Error: Ring is not set.')
            sys.exit()
        self.fs_vis=f_vis_calc(r=self.rs, lb=self.lb, m=self.m) #[c/Rg]
        self.f_vis_set_done=True
    
    def v_rad_set(self):
        if self.f_vis_set_done==False:
            print('Error: viscous frequency is not set.')
            sys.exit()
        self.vs_rad=self.rs*self.fs_vis #[c]

    def epsilon_set(self):
        if self.epsilon_par_set_done==False:
            print('Error: parameter is not set.')
            sys.exit()
        self.eps=epsilon_calc(r=self.rs, stress=self.stress, gamma=self.gamma_eps, r_min=self.r_min_eps)
        self.eps*=2.*np.pi*self.rs*self.wids
        self.eps_tot=np.sum(self.eps)

    # Final operation!
    def out2in(self):
        self.rs=self.rs[::-1]
        self.rs_min=self.rs_min[::-1]
        self.wids=self.wids[::-1]
        self.drs=self.drs[::-1]
        self.fs_vis=self.fs_vis[::-1]
        self.vs_rad=self.vs_rad[::-1]
        self.eps=self.eps[::-1]

# Kepler frequency
def f_kep_calc(r):
    f_kep=(r**(-1.5))/(2.*np.pi) #[c/rg]
    return f_kep

# Viscous frequency (empirical)
def f_vis_calc(r, lb, m):
    f_kep=f_kep_calc(r=r)
    f_vis=lb*(r**(-m))*f_kep #[c/rg]
    return f_vis

# Emissivity
def epsilon_calc(r, stress, gamma, r_min):
    #Stressed
    if stress==True:
        epsilon=r**(-gamma)

    #Stress-free
    elif stress==False:
        epsilon=3.*(1.-np.sqrt(r_min/r))*(r**(-gamma))

    else:
        print('Error')
        sys.exit()

    return epsilon

# ------------------------------------ #
# ----- Mass accretion rate PSD ------ #
# ------------------------------------ #
class FluPro:
    def __init__(self):
        self.mu=1. # fixed without the loss of generarity
        self.sigma=0
        self.psd_wo_prop_done=False
        self.norm_psd=0.
        self.n_data=0
        self.dt=0
        self.f_set_done=False
        self.f_vis_set_done=False

    def f_set(self, fs):
        self.fs=fs #[Hz]
        self.f_set_done=True

    def f_vis_set(self, fs_vis):
        self.fs_vis=fs_vis #[c_rg]
        self.f_vis_set_done=True

    # Calculate PSD without propagation
    def psd_wo_prop(self, c_rg):
        if self.f_set_done==False:
            print('Error: parameter is not set.')
            sys.exit()
        if self.f_vis_set_done==False:
            print('Error: parameter is not set.')
            sys.exit()
        if self.sigma==0:
            print('Error: parameter is not set.')
            sys.exit()
        for i, f_vis in enumerate(self.fs_vis):
            df=f_vis*c_rg # [Hz]
            # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2
            psd_intr=lorentz(f=self.fs, mu=self.mu, sigma=self.sigma, f_c=0, df=df) 

            if i==0:
                self.psds_intr=psd_intr
            else:
                self.psds_intr=np.vstack((self.psds_intr, psd_intr))
        self.psd_wo_prop_done=True

    # Calculate PSD with propagation
    def psd_w_prop(self):
        if self.psd_wo_prop_done==False:
            print('Error: PSD without propagation is not calculated.')
            sys.exit()

        self.norm_psd=2.*self.dt/((self.mu**2)*self.n_data)
        for i in range(len(self.fs_vis)):
            b2s=(1./self.norm_psd)*self.psds_intr[i] # Modulus square of the Fourier transform 
            if i==0:
                #lm2s=self.mdot0*b2s
                lm2s=b2s # Modulus square of the Fourier transform
                psd_prop=self.norm_psd*lm2s # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2
                self.psds_prop=psd_prop
                psd_prop_pre=psd_prop

            else:
                c2s=(1./self.norm_psd)*psd_prop_pre # Modulus square of the Fourier transform
                fs_exte, b2s_exte=ft_mod2_unfold(fs=self.fs, b2s=b2s, mu=self.mu)
                fs_exte, c2s_exte=ft_mod2_unfold(fs=self.fs, b2s=c2s, mu=self.mu)

                lm2s_exte=conv_calc_eff(bs=b2s_exte, cs=c2s_exte)
                fs, lm2s=ft_mod2_fold(fs=fs_exte, bs=lm2s_exte)
                # Modulus square of the Fourier transform
                # |X_{j}|^2=( (|A_{j}|^2) \odot (|B_{j}|^2) ) / N^2 in our definition of Fourier transform
                lm2s/=self.n_data**2 
                psd_prop=self.norm_psd*lm2s # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2

                self.psds_prop=np.vstack((self.psds_prop, psd_prop))
                psd_prop_pre=psd_prop

def lorentz(f, mu, sigma, f_c, df):
    var=sigma**2
    # Normalization is always appropriate: \int _0 ^{\infty} df P(f) = (\sigma/\mu)^2 (fractional variance) 
    # (Ingram & Motta, 2019)
    l=(((sigma/mu)**2)/((np.pi/2.)+np.arctan(f_c/df)))*(df/( ((f-f_c)**2) + (df**2) )) 
    return l

# ---------------------------------------- #
# ----- Mass accretion rate to Flux ------ #
# ---------------------------------------- #
class Mdot2Flux:
    def __init__(self):
        self.e_min=0. # Energy [keV]
        self.e_max=0. # Energy [keV]
        self.n_sc=0 # Number of spectral components
        self.la_eff=0. # Effective area [cm^2]
        self.ws=0.
        self.ws_ref=0.
        self.mu_fl=0.
        self.mu_fl_ref=0.
        self.norm_rep=0.
        self.norm_rep_ref=0.
        self.norm_psd=0.
        self.norm_csd=0.
        self.par_set_done=False
        self.csd_flux_calc_done=False

    def ene_set(self, e_min, e_max):
        self.e_min=e_min
        self.e_max=e_max

    def par_set(self):
        self.fl=np.zeros(self.n_sc)
        self.par_set_done=True

    def la_eff_set(self, es_min, es_max, specresp):
        for i_e, (e_min, e_max) in enumerate(zip(es_min, es_max)):
            if e_min<=self.e and self.e<=e_max:
                self.la_eff=specresp[i_e]
                break
            if i_e==len(es_min)-1:
                print('Error: effective area is not known for E={0:.2f} keV.'.format(self.e))
                sys.exit()

    def flux_set(self, es, specs):
        if self.par_set_done==False:
            print('Error: parameter is not set.')
            sys.exit()
        for i_e, e in enumerate(es):
            if e>self.e:
                for i_sc in range(self.n_sc):
                    self.fl[i_sc]=specs[i_sc][i_e]
                break

    def weight_mu_set(self,\
                      es_area_min,\
                      es_area_max,\
                      las_eff,\
                      es_spec,\
                      spec_disk,\
                      spec_scomp,\
                      spec_hcomp,\
                      spec_sref,\
                      spec_href,\
                      eps_disk,\
                      eps_scomp,\
                      eps_hcomp,\
                      eps_tot_disk,\
                      eps_tot_scomp,\
                      eps_tot_hcomp,\
                      f_dir):

        ws_tot_e,\
        speceff_disk_tot,\
        speceff_sref_tot,\
        speceff_href_tot,\
        speceff_tot,\
        de_tot\
        =weight_spec_sum(e_min=self.e_min,\
                         e_max=self.e_max,\
                         es_area_min=es_area_min,\
                         es_area_max=es_area_max,\
                         las_eff=las_eff,\
                         es_spec=es_spec,\
                         spec_disk=spec_disk,\
                         spec_scomp=spec_scomp,\
                         spec_hcomp=spec_hcomp,\
                         spec_sref=spec_sref,\
                         spec_href=spec_href,\
                         eps_disk=eps_disk,\
                         eps_scomp=eps_scomp,\
                         eps_hcomp=eps_hcomp,\
                         eps_tot_disk=eps_tot_disk,\
                         eps_tot_scomp=eps_tot_scomp,\
                         eps_tot_hcomp=eps_tot_hcomp,\
                         f_dir=f_dir)

        # --- Averaging over energy range of interest --- #
        ws=ws_tot_e/de_tot
        speceff_disk=speceff_disk_tot/de_tot
        speceff_sref=speceff_sref_tot/de_tot
        speceff_href=speceff_href_tot/de_tot
        speceff=speceff_tot/de_tot

        # --- Set quantities necesssary for calculating timing properties --- #
        self.ws=ws
        self.speceff_disk=speceff_disk
        self.speceff_sref=speceff_sref
        self.speceff_href=speceff_href
        self.w_tot=np.sum(self.ws) # Unnecessary? (2021/08/17)
        # <x(E, t)>=(S_{hc}(E)+S_{sc}(E)+S_{d}(E))*A_{eff}(E) 
        # (Absorption should be included in S_{xx}(E))
        # (<x(E, t)> is not affected by f_{dir})
        self.mu_fl=speceff 

    def psd_norm_set(self, dt, n_data):
        # PSD [(sigma/mean)^2/Hz]
        self.norm_psd=2.*dt/((self.mu_fl**2)*n_data) 

    def csd_norm_set(self, dt, n_data):
        # CSD [(sigma/mean)^2/Hz]
        self.norm_csd=2.*dt/((self.mu_fl_ref*self.mu_fl)*n_data) 

    #def norm_rep_set(self, f_rep, dt0, w_flow_tot):
    #    #C(E) in the impulse response
    #    self.norm_rep=(f_rep*self.speceff_disk+self.speceff_sref+self.speceff_href)/(w_flow_tot*dt0) 

    def norm_rep_set(self, dt0, w_flow_tot):
        #C(E) in the impulse response
        ref=self.speceff_sref+self.speceff_mref+self.speceff_href
        self.norm_rep=ref/(w_flow_tot*dt0) 

    #def norm_rep_ref_set(self, f_rep, dt0, w_flow_tot):
    #    #C(E) in the impulse response
    #    self.norm_rep_ref=(f_rep*self.speceff_disk_ref+self.speceff_sref_ref+self.speceff_href_ref)/(w_flow_tot*dt0) 

    def norm_rep_ref_set(self, dt0, w_flow_tot):
        #C(E) in the impulse response
        ref=self.speceff_sref_ref+self.speceff_mref_ref+self.speceff_href_ref
        self.norm_rep=ref/(w_flow_tot*dt0) 

    # Reprocessing is included.
    def psd_flux_rep_calc(self,\
                          fs,\
                          n_r,\
                          ws_flow,\
                          lm2s,\
                          fs_vis,\
                          dr_r,\
                          t0,\
                          dt0,\
                          rg_c):
        if self.norm_psd==0:
            print('Error: normalization of PSD is not set.')
            sys.exit()

        #|X(f)|^2, Direct component
        lf2s_dir=lf2_calc(fs=fs,\
                          n_r=n_r,\
                          ws=self.ws,\
                          lm2s=lm2s,\
                          fs_vis=fs_vis,\
                          dr_r=dr_r,\
                          rg_c=rg_c) 

        #|X(f)|^2, (one term of) Reprocessed component
        lf2s_rep=lf2_calc(fs=fs,\
                          n_r=n_r,\
                          ws=ws_flow,\
                          lm2s=lm2s,\
                          fs_vis=fs_vis,\
                          dr_r=dr_r,\
                          rg_c=rg_c)*\
                 lh2_rep_calc(f=fs, norm=self.norm_rep, dt0=dt0)

        #|(X(f))^{*}Y(f)|^2, (one term of) Reprocessed component (cross term)
        lflfs=lflf_calc(fs=fs,\
                        n_r=n_r,\
                        ws_ref=self.ws,\
                        ws_coi=ws_flow,\
                        lm2s=lm2s,\
                        fs_vis=fs_vis,\
                        dr_r=dr_r,\
                        rg_c=rg_c)*\
              lh_rep_calc(f=fs, norm=self.norm_rep, t0=t0, dt0=dt0)
        lflfs=2.*lflfs.real

        # Total
        lf2s=lf2s_dir+lf2s_rep+lflfs
        # Normalize
        self.psd_fl=self.norm_psd*lf2s

    # Reprocessing is included.
    def csd_flux_rep_calc(self,\
                          fs,\
                          n_r,\
                          ws_rep,\
                          lm2s,\
                          fs_vis,\
                          dr_r,\
                          t0,\
                          dt0,\
                          rg_c):
        if self.norm_csd==0:
            print('Error: normalization of CSD is not set.')
            sys.exit()
        self.fs=fs

        #|(X(f))^{*}Y(f)|^2, Direct vs Direct
        csd_dirdir=lflf_calc(fs=self.fs,\
                             n_r=n_r,\
                             ws_ref=self.ws_ref,\
                             ws_coi=self.ws,\
                             lm2s=lm2s,\
                             fs_vis=fs_vis,\
                             dr_r=dr_r,\
                             rg_c=rg_c)

        #|(X(f))^{*}Y(f)|^2, Reprocessed vs Reprocessed
        csd_reprep=lf2_calc(fs=fs,\
                            n_r=n_r,\
                            ws=ws_rep,\
                            lm2s=lm2s,\
                            fs_vis=fs_vis,\
                            dr_r=dr_r,\
                            rg_c=rg_c)*\
                   (lh_rep_calc(f=self.fs, norm=self.norm_rep_ref, t0=t0, dt0=dt0).conjugate())*\
                    lh_rep_calc(f=self.fs, norm=self.norm_rep, t0=t0, dt0=dt0)

        #|(X(f))^{*}Y(f)|^2, Direct vs Reprocessed
        csd_dirrep=lflf_calc(fs=self.fs,\
                             n_r=n_r,\
                             ws_ref=self.ws_ref,\
                             ws_coi=ws_rep,\
                             lm2s=lm2s,\
                             fs_vis=fs_vis,\
                             dr_r=dr_r,\
                             rg_c=rg_c)*\
                   lh_rep_calc(f=self.fs, norm=self.norm_rep, t0=t0, dt0=dt0)

        #|(X(f))^{*}Y(f)|^2, Reprocessed vs Direct
        csd_repdir=lflf_calc(fs=self.fs,\
                             n_r=n_r,\
                             ws_ref=ws_rep,\
                             ws_coi=self.ws,\
                             lm2s=lm2s,\
                             fs_vis=fs_vis,\
                             dr_r=dr_r,\
                             rg_c=rg_c)*\
                   (lh_rep_calc(f=self.fs, norm=self.norm_rep_ref, t0=t0, dt0=dt0).conjugate())

        # Total
        self.csd_fl=csd_dirdir+\
                    csd_reprep+\
                    csd_dirrep+\
                    csd_repdir
        # Normalize
        self.csd_fl=self.norm_csd*self.csd_fl
        self.csd_flux_calc_done=True

    def lagf_calc(self):
        if self.csd_flux_calc_done==False:
            print('Error: cross spectrum is not set.')
            sys.exit()

        # ----- Phase lag ----- #
        for i_ele, csd_fl_ele in enumerate(self.csd_fl):
            phi_ele=cmath.phase(csd_fl_ele)
            if i_ele==0:
                self.phi=np.array(phi_ele)
            else:
                self.phi=np.append(self.phi, phi_ele)

        # ----- Time lag ----- #
        self.tau=self.phi/(2.*np.pi*self.fs) #Positive lag = hard lag

# ---------------------------------- #
# ---------- Reprocessing ---------- #
# ---------------------------------- #
# Energy dependent normalization of the impulse response for reprocessing
def norm_h_dc_calc(frac_dc, spec_hcomp, spec_scomp, dt0):
    norm=frac_dc*(spec_hcomp+spec_scomp)/dt0
    return norm

# Impulse response for reprocessing
def h_rep_calc(t, norm, t0, dt0):
    h=norm*(t0-dt0/2.<t)*(t<t0+dt0/2.)
    return h

# Transfer function for reprocessing
def lh_rep_calc(f, norm, t0, dt0):
    omega=2.*np.pi*f
    #lh=norm*dt0*np.sinc(omega*dt0/2.)*np.exp(-1j*omega*t0) # Wrong!
    lh=norm*dt0*np.sinc(omega*dt0/(2.*np.pi))*np.exp(-1j*omega*t0) # Correct!
    return lh

# Modulus square of transfer function for reprocessing
def lh2_rep_calc(f, norm, dt0):
    omega=2.*np.pi*f
    #lh2=(norm*dt0*np.sinc(omega*dt0/2.))**2 # Wrong!
    lh2=(norm*dt0*np.sinc(omega*dt0/(2.*np.pi)))**2 # Correct!
    return lh2

# ------------------------------------------------ #
# ---------- Calculation of PSD and CSD ---------- #
# ------------------------------------------------ #
def prop_time_calc(i_start, i_end, ts_vis, dr_r):
    t_prop=0.
    for i in range(i_start, i_end):
        t_prop+=ts_vis[i]
    t_prop*=dr_r
    return t_prop

def lf2_calc(fs,\
             n_r,\
             ws,\
             lm2s,\
             fs_vis,\
             dr_r,\
             rg_c):
    ts_vis=1./fs_vis #[Rg/c]
    tot=0
    for i_r in range(n_r):
        tot+=(ws[i_r]**2)*lm2s[i_r]

        if i_r==0:
            continue

        # Cross term
        tot_c=0
        for i_ro in range(i_r):
            # Propagation time
            t_prop=prop_time_calc(i_start=i_ro, i_end=i_r, ts_vis=ts_vis, dr_r=dr_r) #[Rg/c]
            t_prop*=rg_c #[s]
            # Cross term
            tot_c+=ws[i_ro]*ws[i_r]*np.cos(2.*np.pi*fs*t_prop)*lm2s[i_ro]

        tot+=2.*tot_c

    return tot

def lflf_calc(fs,\
              n_r,\
              ws_ref,\
              ws_coi,\
              lm2s,\
              fs_vis,\
              dr_r,\
              rg_c):
    ts_vis=1./fs_vis #[Rg/c]
    tot=0
    for i_r in range(n_r):
        tot+=ws_ref[i_r]*ws_coi[i_r]*lm2s[i_r]

        if i_r==0:
            continue

        # Cross term
        tot_c=0
        for i_ro in range(i_r):
            # Propagation time
            t_prop=prop_time_calc(i_start=i_ro, i_end=i_r, ts_vis=ts_vis, dr_r=dr_r) #[Rg/c]
            t_prop*=rg_c #[s]
            # Cross term
            # Ingram & van der Klis
            #tot_c+=(ws_ref[i_ro]*ws_coi[i_r]*np.exp(1j*2.*np.pi*fs*t_prop) + \
            #        ws_ref[i_r]*ws_coi[i_ro]*np.exp(-1j*2.*np.pi*fs*t_prop))*lm2s[i_ro]

            # Mofification due to the difference of the Fourier transform
            tot_c+=(ws_ref[i_ro]*ws_coi[i_r]*np.exp(-1j*2.*np.pi*fs*t_prop) + \
                    ws_ref[i_r]*ws_coi[i_ro]*np.exp(1j*2.*np.pi*fs*t_prop))*lm2s[i_ro]

        tot=tot+tot_c

    return tot

def print_ring_info(name, xs, digit):
    print('--------------------------------------------------')
    print('{0}: '.format(name), end='')
    n_x=len(xs)
    for i_x, x in enumerate(xs):
        if i_x==n_x-1:
            print_digit_end(x=x, digit=digit, end=True)
        else:
            print_digit_end(x=x, digit=digit, end=False)

def print_digit_end(x, digit, end):
    if digit==1:
        if end==False:
            print('{:.1f}, '.format(x), end='')
        else:
            print('{:.1f}'.format(x))
    if digit==2:
        if end==False:
            print('{:.2f}, '.format(x), end='')
        else:
            print('{:.2f}'.format(x))
    if digit==3:
        if end==False:
            print('{:.3f}, '.format(x), end='')
        else:
            print('{:.3f}'.format(x))
    if digit==4:
        if end==False:
            print('{:.4f}, '.format(x), end='')
        else:
            print('{:.4f}'.format(x))
    if digit==5:
        if end==False:
            print('{:.5f}, '.format(x), end='')
        else:
            print('{:.5f}'.format(x))

