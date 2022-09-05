import sys
import cmath
import xspec
from fourier_h import *

# -------------------------- #
# ----- Read parameter ----- #
# -------------------------- #
class SetParameter:
    def __init__(self):
        self.set_inpar_done=False

    # --- Model: preflow --- #
    def preflow_set_inpar(self, pars, es):
        self.mass        =pars[0]  # BH mass [solar mass]
        self.r_in        =pars[1]  # Inner radius of hard Compton [Rg]
        self.dr_hcomp    =pars[2]  # Width of hard Compton [Rg]
        self.dr_scomp    =pars[3]  # Width of soft Compton [Rg]
        self.dr_disk     =pars[4]  # Width of variable disk [Rg]
        self.n_ring      =pars[5]  # Number of rings splitted [Rg]
        self.cf_var_d    =pars[6]  # Fractional variability at variable disk [-]
        self.cf_var_s    =pars[7]  # Fractional variability at soft Compton [-]
        self.cf_var_h    =pars[8]  # Fractional variability at hard Compton [-]
        self.cb_d        =pars[9]  # B_{gen}_{disk} [-]
        self.m_d         =pars[10] # m_{gen}_{disk} [-]
        self.cbp_d       =pars[11] # B_{gen}_{disk} [-]
        self.mp_d        =pars[12] # m_{gen}_{disk} [-]
        self.cb_f        =pars[13] # B_{prop}_{flow} [-]
        self.m_f         =pars[14] # m_{prop}_{flow} [-]
        self.cbp_f       =pars[15] # B_{prop}_{flow} [-]
        self.mp_f        =pars[16] # m_{prop}_{flow} [-]
        self.cs          =pars[17] # Damping factor [-]
        self.gamma_disk  =pars[18] # Radial index of emissivity [-]
        self.gamma_flow  =pars[19] # Radial index of emissivity [-]
        self.stress      =pars[20] # 1: stressed, 2: stress-free in emissivity
        self.e_min       =pars[21] # Lower bound of energy band [keV] (unused)
        self.e_max       =pars[22] # Upper bound of energy band [keV] (unused)
        self.cs_d        =pars[23] # Disk spectrum [-]
        self.cs_s        =pars[24] # Soft Compton spectrum [-]
        self.cs_h        =pars[25] # Hard Compton spectrum [-]
        self.cs_sr       =pars[26] # Soft reflection spectrum [-]
        self.cs_hr       =pars[27] # Hard Compton spectrum [-]
        self.csp_s       =pars[28] # Soft Compton sensitive to seed photon variability [-]
        self.csp_h       =pars[29] # Hard Compton sensitive to seed photon variability [-]
        self.csp_sr      =pars[30] # Soft reflection sensitive to seed photon variability [-]
        self.csp_hr      =pars[31] # Hard reflection sensitive to seed photon variability [-]
        self.e_minr      =pars[32] # Lower bound of reference band [keV] (unused)
        self.e_maxr      =pars[33] # Upper bound of reference band [keV] (unused)
        self.cs_d_r      =pars[34] # Disk spectrum [-]
        self.cs_s_r      =pars[35] # Soft Compton spectrum [-]
        self.cs_h_r      =pars[36] # Hard Compton spectrum [-]
        self.cs_sr_r     =pars[37] # Soft reflection spectrum [-]
        self.cs_hr_r     =pars[38] # Hard reflection spectrum [-]
        self.csp_s_r     =pars[39] # Soft Compton spectrum [-]
        self.csp_h_r     =pars[40] # Hard Compton spectrum [-]
        self.csp_sr_r    =pars[41] # Soft reflection spectrum [-]
        self.csp_hr_r    =pars[42] # Hard reflection spectrum [-]
        self.tr_s        =pars[43] # Rising time of impulse response for soft reflection [sec]
        self.dt0_s       =pars[44] # Time width of impulse response for soft reflection [sec]
        self.tr_h        =pars[45] # Rising time of impulse response for hard reflection [sec]
        self.dt0_h       =pars[46] # Time width of impulse response for hard reflection  [sec]
        self.quant       =pars[47]
            # 1: power spectrum 
            # 2: real part of cross spectrum
            # 3: imaginary part of cross spectrum
            # 4: absolute value of cross spectrum
            # 5: phase lag (Positive lag means reference band lagging behind energy band.)
            # 6: time lag  (Positive lag means reference band lagging behind energy band.)
        self.invert      =pars[48] # 1: Normal,  2: Im[C(f)], phase lag, and time lag are multiplied by -1.
        self.display     =pars[49] # 1: display, 2: not display

        # PREFLOW model is a timing model!
        # Energy in XSPEC corresponds to Fourier frequency in preflow.
        self.fs_data=es

        # Hidden parameters 
        self.r_min       =self.r_in # Radial index of emissivity [-]
        self.tr_d         =0.    # Time width of reflection impulse response (variable disk) [sec]
        self.dt0_d        =1.e-2 # Time width of reflection impulse response (variable disk) [sec]

        ### Geometry ###
        self.r_sh        =self.r_in+self.dr_hcomp  # Inner radius of soft Compton [Rg]
        self.r_ds        =self.r_sh+self.dr_scomp  # Inner radius of disk [Rg]
        self.r_out       =self.r_ds+self.dr_disk   # Outer radius of disk [Rg]

        # Impulse response of reflection
        self.t0_d=self.tr_d+(self.dt0_d/2.)
        self.t0_s=self.tr_s+(self.dt0_s/2.)
        self.t0_h=self.tr_h+(self.dt0_h/2.)

        # For the case that the disk is not variable
        if self.cf_var_d==0.:
            self.r_out=self.r_ds

        self.set_inpar_done=True

    # --- Model: preflows --- #
    def preflows_set_inpar(self, pars, es):
        ###########################
        ### Spectral parameters ###
        ###########################
        # Disk (diskbb is employed)
        self.temp_bb_d    =pars[0]  # Blackbody temperature [keV]
        self.norm_d       =pars[1]  # Normalization of disk [-]
        # Soft Compton (cutoffpl is employed)
        self.gamma_s      =pars[2]  # Spectral index of soft Compton [-]
        self.ecut_s       =pars[3]  # Electron temperature of soft Compton [keV]
        self.norm_s       =pars[4]  # Normalization of soft Compton [-]
        # Hard Compton (cutoffpl is employed)
        self.gamma_h      =pars[5]  # Spectral index of hard Compton [-]
        self.ecut_h       =pars[6]  # Electron temperature of hard Compton [keV]
        self.norm_h       =pars[7]  # Normalization of hard Compton [-]
        # Reflection (common) (relxill is employed)
        self.incl         =pars[8]
        self.a            =pars[9]
        self.cafe         =pars[10]
        # Soft reflection
        self.r_in_sr      =pars[11]
        self.r_out_sr     =pars[12]
        self.index_sr     =pars[13]
        self.logxi_sr     =pars[14]
        self.cf_sr        =pars[15] # This is normalization, not refl_frac!
        # Hard reflection
        self.r_in_hr      =pars[16]
        self.r_out_hr     =pars[17]
        self.index_hr     =pars[18]
        self.logxi_hr     =pars[19]
        self.cf_hr        =pars[20] # This is normalization, not refl_frac! 

        ##############################
        ### Variability parameters ###
        ##############################
        self.mass         =pars[21]  # BH mass [solar mass]
        self.r_in         =pars[22]  # Inner radius of hard Compton (inner radius of hot flow) [Rg]
        self.dr_hcomp     =pars[23]  # Inner radius of hard Compton (inner radius of hot flow) [Rg]
        self.dr_scomp     =pars[24]  # Transition radius between soft Compton and variable disk [Rg]
        self.dr_disk      =pars[25]  # Transition radius between soft Compton and variable disk (outer radius of hot flow) [Rg]
        self.n_ring       =pars[26]  # Outer radius of variable disk [Rg]
        self.cf_var_d     =pars[27]  # Fractional variability at variable disk [-]
        self.cf_var_s     =pars[28]  # Fractional variability at soft Compton [-]
        self.cf_var_h     =pars[29]  # Fractional variability at hard Compton [-]
        self.cb_d         =pars[30]  # B_{disk} [-]
        self.m_d          =pars[31]  # m_{disk} [-]
        self.cbp_d        =pars[32]  # Bprop_{disk} [-]
        self.mp_d         =pars[33]  # mprop_{disk} [-]
        self.cb_f         =pars[34]  # B_{flow} [-]
        self.m_f          =pars[35]  # m_{flow} [-]
        self.cbp_f        =pars[36]  # Bprop_{flow} [-]
        self.mp_f         =pars[37]  # mprop_{flow} [-]
        self.cs           =pars[38]  # Smoothing factor [-]
        self.gamma_disk   =pars[39]  # Radial index of emissivity [-]
        self.gamma_flow   =pars[40]  # Radial index of emissivity [-]
        self.stress       =pars[41]  # 1: stressed, 2: stress-free in emissivity
        self.e_min        =pars[42]  # Lower bound of energy band [keV] (unused)
        self.e_max        =pars[43]  # Upper bound of energy band [keV] (unused)
        self.e_minr       =pars[44]  # Lower bound of reference band [keV] (unused)
        self.e_maxr       =pars[45]  # Upper bound of reference band [keV] (unused)
        self.eta0_d       =pars[46]
        self.eta1_d       =pars[47]
        self.eta0_s       =pars[48]
        self.eta1_s       =pars[49]
        self.eta0_h       =pars[50]
        self.eta1_h       =pars[51]
        self.etap0_s      =pars[52] # Sensitivity parameter for seed photon variability
        self.etap1_s      =pars[53] # Sensitivity parameter for seed photon variability
        self.etap0_h      =pars[54] # Sensitivity parameter for seed photon variability
        self.etap1_h      =pars[55] # Sensitivity parameter for seed photon variability
        self.tr_s         =pars[56] # Start time of reflection impulse response [sec]
        self.dt0_s        =pars[57] # Time width of reflection impulse response [sec]
        self.tr_h         =pars[58] # Start time of reflection impulse response [sec]
        self.dt0_h        =pars[59] # Time width of reflection impulse response [sec]
        self.quant        =pars[60]
            # 0: energy spectrum
            # 1: power spectrum 
            # 2: real part of cross spectrum
            # 3: imaginary part of cross spectrum
            # 4: absolute value of cross spectrum
            # 5: phase lag (Positive lag means reference band lagging behind energy band.)
            # 6: time lag  (Positive lag means reference band lagging behind energy band.)
        self.invert       =pars[61] # 1: Normal,  2: Im[C(f)], phase lag, and time lag are multiplied by -1.
        self.display      =pars[62] # 1: display, 2: not display

        # PREFLOW model is a spectral-timing model!
        # Energy in XSPEC corresponds to Fourier frequency in preflow.
        self.es=np.array(es)
        self.fs_data=es

        # Hidden parameters 
        self.r_min        =self.r_in  # Minumum radius of emissivity [-]
        self.tr_d         =0.    # Time width of reflection impulse response (variable disk) [sec]
        self.dt0_d        =1.e-2 # Time width of reflection impulse response (variable disk) [sec]
        self.eta0_sr      =self.eta0_s
        self.eta1_sr      =self.eta1_s
        self.eta0_hr      =self.eta0_h
        self.eta1_hr      =self.eta1_h
        self.etap0_sr     =self.etap0_s
        self.etap1_sr     =self.etap1_s
        self.etap0_hr     =self.etap0_h
        self.etap1_hr     =self.etap1_h

        ### Geometry ###
        self.r_sh        =self.r_in+self.dr_hcomp  # Inner radius of soft Compton [Rg]
        self.r_ds        =self.r_sh+self.dr_scomp  # Inner radius of disk [Rg]
        self.r_out       =self.r_ds+self.dr_disk   # Outer radius of disk [Rg]

        # Impulse response of reflection
        self.t0_d=self.tr_d+(self.dt0_d/2.)
        self.t0_s=self.tr_s+(self.dt0_s/2.)
        self.t0_h=self.tr_h+(self.dt0_h/2.)

        # For the case that the disk is not variable
        if self.cf_var_d==0.:
            self.r_out=self.r_ds

        self.set_inpar_done=True

    # --- Model: preflowscp --- #
    def preflowscp_set_inpar(self, pars, es):
        ###########################
        ### Spectral parameters ###
        ###########################
        # Disk (diskbb is employed)
        self.temp_bb_d    =pars[0]  # Blackbody temperature [keV]
        self.norm_d       =pars[1]  # Normalization of disk [-]
        # Soft Compton (nthcomp is employed)
        self.gamma_s      =pars[2]  # Spectral index of soft Compton [-]
        self.temp_bb_s    =pars[3]  # Seed photon temperature of soft Compton [keV]
        self.temp_e_s     =pars[4]  # Electron temperature of soft Compton [keV]
        self.norm_s       =pars[5]  # Normalization of soft Compton [-]
        # Hard Compton (nthcomp is employed)
        self.gamma_h      =pars[6]  # Spectral index of hard Compton [-]
        self.temp_bb_h    =pars[7]  # Seed photon temperature of hard Compton [keV]
        self.temp_e_h     =pars[8]  # Electron temperature of hard Compton [keV]
        self.norm_h       =pars[9]  # Normalization of hard Compton [-]
        # Reflection (common) (relxillCp is employed)
        self.incl         =pars[10]
        self.a            =pars[11]
        self.cafe         =pars[12]
        # Soft reflection
        self.r_in_sr      =pars[13]
        self.r_out_sr     =pars[14]
        self.index_sr     =pars[15]
        self.logxi_sr     =pars[16]
        self.logcn_sr     =pars[17]
        self.cf_sr        =pars[18] # This is normalization, not refl_frac!
        # Hard reflection
        self.r_in_hr      =pars[19]
        self.r_out_hr     =pars[20]
        self.index_hr     =pars[21]
        self.logxi_hr     =pars[22]
        self.logcn_hr     =pars[23]
        self.cf_hr        =pars[24] # This is normalization, not refl_frac! 

        ##############################
        ### Variability parameters ###
        ##############################
        self.mass         =pars[25]  # BH mass [solar mass]
        self.r_in         =pars[26]  # Inner radius of hard Compton (inner radius of hot flow) [Rg]
        self.dr_hcomp     =pars[27]  # Inner radius of hard Compton (inner radius of hot flow) [Rg]
        self.dr_scomp     =pars[28]  # Transition radius between soft Compton and variable disk [Rg]
        self.dr_disk      =pars[29]  # Transition radius between soft Compton and variable disk (outer radius of hot flow) [Rg]
        self.n_ring       =pars[30]  # Outer radius of variable disk [Rg]
        self.cf_var_d     =pars[31]  # Fractional variability at variable disk [-]
        self.cf_var_s     =pars[32]  # Fractional variability at soft Compton [-]
        self.cf_var_h     =pars[33]  # Fractional variability at hard Compton [-]
        self.cb_d         =pars[34]  # B_{disk} [-]
        self.m_d          =pars[35]  # m_{disk} [-]
        self.cbp_d        =pars[36]  # Bprop_{disk} [-]
        self.mp_d         =pars[37]  # mprop_{disk} [-]
        self.cb_f         =pars[38]  # B_{flow} [-]
        self.m_f          =pars[39]  # m_{flow} [-]
        self.cbp_f        =pars[40]  # Bprop_{flow} [-]
        self.mp_f         =pars[41]  # mprop_{flow} [-]
        self.cs           =pars[42]  # Smoothing factor [-]
        self.gamma_disk   =pars[43]  # Radial index of emissivity [-]
        self.gamma_flow   =pars[44]  # Radial index of emissivity [-]
        self.stress       =pars[45]  # 1: stressed, 2: stress-free in emissivity
        self.e_min        =pars[46]  # Lower bound of energy band [keV] (unused)
        self.e_max        =pars[47]  # Upper bound of energy band [keV] (unused)
        self.e_minr       =pars[48]  # Lower bound of reference band [keV] (unused)
        self.e_maxr       =pars[49]  # Upper bound of reference band [keV] (unused)
        self.eta0_d       =pars[50]
        self.eta1_d       =pars[51]
        self.eta0_s       =pars[52]
        self.eta1_s       =pars[53]
        self.eta0_h       =pars[54]
        self.eta1_h       =pars[55]
        self.etap0_s      =pars[56] # Sensitivity parameter for seed photon variability
        self.etap1_s      =pars[57] # Sensitivity parameter for seed photon variability
        self.etap0_h      =pars[58] # Sensitivity parameter for seed photon variability
        self.etap1_h      =pars[59] # Sensitivity parameter for seed photon variability
        self.tr_s         =pars[60] # Start time of reflection impulse response [sec]
        self.dt0_s        =pars[61] # Time width of reflection impulse response [sec]
        self.tr_h         =pars[62] # Start time of reflection impulse response [sec]
        self.dt0_h        =pars[63] # Time width of reflection impulse response [sec]
        self.quant        =pars[64]
            # 0: energy spectrum
            # 1: power spectrum 
            # 2: real part of cross spectrum
            # 3: imaginary part of cross spectrum
            # 4: absolute value of cross spectrum
            # 5: phase lag (Positive lag means reference band lagging behind energy band.)
            # 6: time lag  (Positive lag means reference band lagging behind energy band.)
        self.invert       =pars[65] # 1: Normal,  2: Im[C(f)], phase lag, and time lag are multiplied by -1.
        self.display      =pars[66] # 1: display, 2: not display

        # PREFLOW model is a spectral-timing model!
        # Energy in XSPEC corresponds to Fourier frequency in preflow.
        self.es=np.array(es)
        self.fs_data=es

        # Hidden parameters 
        self.r_min        =self.r_in  # Minumum radius of emissivity [-]
        self.tr_d         =0.    # Time width of reflection impulse response (variable disk) [sec]
        self.dt0_d        =1.e-2 # Time width of reflection impulse response (variable disk) [sec]
        self.eta0_sr      =self.eta0_s
        self.eta1_sr      =self.eta1_s
        self.eta0_hr      =self.eta0_h
        self.eta1_hr      =self.eta1_h
        self.etap0_sr     =self.etap0_s
        self.etap1_sr     =self.etap1_s
        self.etap0_hr     =self.etap0_h
        self.etap1_hr     =self.etap1_h

        ### Geometry ###
        self.r_sh        =self.r_in+self.dr_hcomp  # Inner radius of soft Compton [Rg]
        self.r_ds        =self.r_sh+self.dr_scomp  # Inner radius of disk [Rg]
        self.r_out       =self.r_ds+self.dr_disk   # Outer radius of disk [Rg]

        # Impulse response of reflection
        self.t0_d=self.tr_d+(self.dt0_d/2.)
        self.t0_s=self.tr_s+(self.dt0_s/2.)
        self.t0_h=self.tr_h+(self.dt0_h/2.)

        # The case the disk is not variable
        if self.cf_var_d==0.:
            self.r_out=self.r_ds

        self.set_inpar_done=True

# ------------------------------------- #
# ----- Analytical spectral model ----- #
# ------------------------------------- #
class SpectralModel:
    def __init__(self):
        pass    

    def disk_spec(self, e, temp, norm):
        ratio=(e/temp)*(e/temp<100.)+100.*(e/temp>=100.) # Upper limit to avoid overflow
        cn=norm*(e**2)/(np.exp(ratio)-1.)
        return cn

    def compton_spec(self, e, gamma, ecut, norm):
        cn=norm*(e**(-gamma))*np.exp(-e/ecut)
        return cn

    def reflection_spec(self, e, gamma, ecut, eline, deline, frac_line, norm):
        cn=norm*( (e**(-gamma))*np.exp(-e/ecut) + frac_line*np.exp(-((e-eline)**2)/(2.*(deline**2))) )
        return cn

# -------------------------------- #
# ----- XSPEC spectral model ----- #
# -------------------------------- #
class SpectralModelXspec:
    def __init__(self):
        pass    

    def tbabs_spec(self, es, nh):
        params=[nh]
        fluxes=[]
        es=es.tolist()
        xspec.callModelFunction(modelName='TBabs', energies=es, params=params, flux=fluxes)
        fluxes=np.array(fluxes)
        return fluxes

    def diskbb_spec(self, es, temp, norm):
        params=[temp]
        fluxes=[]
        es=es.tolist()
        xspec.callModelFunction(modelName='diskbb', energies=es, params=params, flux=fluxes)
        fluxes=np.array(fluxes)
        fluxes*=norm
        return fluxes

    def cutoffpl_spec(self, es, gamma, ecut, norm):
        params=[gamma, ecut]
        fluxes=[]
        es=es.tolist()
        xspec.callModelFunction(modelName='cutoffpl', energies=es, params=params, flux=fluxes)
        fluxes=np.array(fluxes)
        fluxes*=norm
        return fluxes

    # Hidden parameters: inp_type=0, Redshift=0.
    def nthcomp_spec(self, es, gamma, ktbb, kte, norm):
        inp_type=0
        redshift=0
        params=[gamma, kte, ktbb, inp_type, redshift]
        fluxes=[]
        es=es.tolist()
        xspec.callModelFunction(modelName='nthComp', energies=es, params=params, flux=fluxes)
        fluxes=np.array(fluxes)
        fluxes*=norm
        return fluxes

    # Hidden parameters: z=0, refl_frac=-1
    def relxill_spec(self, es, incl, a, r_in, r_out, index, gamma, logxi, cafe, kte, norm):
        if norm==0:
            fluxes=0.*np.ones(len(es)-1)
        else:
            z=0
            r_br=(r_in+r_out)/2.
            # refl_frac (Fref) must be -1 so that relxillCp calculates only reflected emission.
            params=[index, index, r_br, a, incl, r_in, r_out, z, gamma, logxi, cafe, kte, -1.]
            fluxes=[]
            es=es.tolist()
            xspec.callModelFunction(modelName='relxill', energies=es, params=params, flux=fluxes)
            fluxes=np.array(fluxes)
            fluxes*=norm
        return fluxes

    # Hidden parameters: z=0, refl_frac=-1
    def relxillcp_spec(self, es, incl, a, r_in, r_out, index, gamma, logxi, logcn, cafe, kte, norm):
        if norm==0:
            fluxes=0.*np.ones(len(es)-1)
        else:
            z=0
            r_br=(r_in+r_out)/2.
            # refl_frac (Fref) must be -1 so that relxillCp calculates only reflected emission.
            params=[incl, a, r_in, r_out, r_br, index, index, z, gamma, logxi, logcn, cafe, kte, -1.]
            fluxes=[]
            es=es.tolist()
            xspec.callModelFunction(modelName='relxillCp', energies=es, params=params, flux=fluxes)
            fluxes=np.array(fluxes)
            fluxes*=norm
        return fluxes

# ------------------------------ #
# ----- Calculate spectrum ----- #
# ------------------------------ #
class EnergySpectrum:
    def __init__(self):
        pass    

    def preflows_calc_spectra(self, es, pars):
        ktbbd =pars[0] 
        normd =pars[1] 
        gammas=pars[2] 
        ecuts =pars[3] 
        norms =pars[4] 
        gammah=pars[5] 
        ecuth =pars[6] 
        normh =pars[7] 
        incl  =pars[8]
        a     =pars[9]
        afe   =pars[10]
        rins  =pars[11]
        routs =pars[12]
        indexs=pars[13]
        logxis=pars[14]
        normsr=pars[15]
        rinh  =pars[16]
        routh =pars[17]
        indexh=pars[18]
        logxih=pars[19]
        normhr=pars[20]

        specx=SpectralModelXspec()
        fluxes_d =specx.diskbb_spec   (es=es, temp=ktbbd, norm=normd)
        fluxes_s =specx.cutoffpl_spec (es=es, gamma=gammas, ecut=ecuts, norm=norms)
        fluxes_h =specx.cutoffpl_spec (es=es, gamma=gammah, ecut=ecuth, norm=normh)
        fluxes_sr=specx.relxill_spec(\
            es=es,        incl=incl,    a=a,\
            r_in=rins,    r_out=routs,  index=indexs,\
            gamma=gammas, logxi=logxis,\
            cafe=afe,     kte=ecuts,    norm=normsr)
        fluxes_hr=specx.relxill_spec(\
            es=es,        incl=incl,    a=a,\
            r_in=rinh,    r_out=routh,  index=indexh,\
            gamma=gammah, logxi=logxih,\
            cafe=afe,     kte=ecuth,    norm=normhr)

        fluxes=fluxes_d+fluxes_s+fluxes_h+fluxes_sr+fluxes_hr

        return fluxes

    def preflowscp_calc_spectra(self, es, pars):
        ktbbd =pars[0] 
        normd =pars[1] 
        gammas=pars[2] 
        ktbbs =pars[3] 
        ktes  =pars[4] 
        norms =pars[5] 
        gammah=pars[6] 
        ktbbh =pars[7] 
        kteh  =pars[8] 
        normh =pars[9] 
        incl  =pars[10]
        a     =pars[11]
        afe   =pars[12]
        rins  =pars[13]
        routs =pars[14]
        indexs=pars[15]
        logxis=pars[16]
        logns =pars[17]
        normsr=pars[18]
        rinh  =pars[19]
        routh =pars[20]
        indexh=pars[21]
        logxih=pars[22]
        lognh =pars[23]
        normhr=pars[24]

        specx=SpectralModelXspec()
        fluxes_d =specx.diskbb_spec   (es=es, temp=ktbbd, norm=normd)
        fluxes_s =specx.nthcomp_spec  (\
            es=es, gamma=gammas, ktbb=ktbbs, kte=ktes, norm=norms)
        fluxes_h =specx.nthcomp_spec  (\
            es=es, gamma=gammah, ktbb=ktbbh, kte=ktes, norm=normh)
        fluxes_sr=specx.relxillcp_spec(\
            es=es,        incl=incl,    a=a,\
            r_in=rins,    r_out=routs,  index=indexs,\
            gamma=gammas, logxi=logxis, logcn=logns,\
            cafe=afe,     kte=ktes,     norm=normsr)
        fluxes_hr=specx.relxillcp_spec(\
            es=es,        incl=incl,    a=a,\
            r_in=rinh,    r_out=routh,  index=indexh,\
            gamma=gammah, logxi=logxih, logcn=lognh,\
            cafe=afe,     kte=kteh,     norm=normhr)

        fluxes=fluxes_d+fluxes_s+fluxes_h+fluxes_sr+fluxes_hr

        return fluxes

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
        # dt = 2^{-n}
        index=-1
        while True:
            index+=1
            dt=2**(-index)
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
        self.fs_prop=np.zeros(self.n_ring)
        self.vs_prop=np.zeros(self.n_ring)
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
        # Change the actual outer radius rout (2022/01/28)
        # This is for setting weight appropriately.
        self.r_out_act=self.rs_min[-1]*(1.+self.k)

    # Center position of ring
    def centerr_calc(self):
        if self.set==False:
            print('Error: parameter is not set.')
            sys.exit()
        for i in range(self.n_ring):
            if i==self.n_ring-1:
                # 2022/01/28
                #self.rs[i]=(self.rs_min[i]+self.r_out)/2.
                self.rs[i]=(self.rs_min[i]+self.r_out_act)/2.
            else:
                self.rs[i]=(self.rs_min[i]+self.rs_min[i+1])/2.

    # Width of ring
    def width_calc(self):
        if self.set==False:
            print('Error: parameter is not set.')
            sys.exit()
        for i in range(self.n_ring):
            if i==self.n_ring-1:
                # 2022/01/28
                #self.wids[i]=self.r_out-self.rs_min[i]
                self.wids[i]=self.r_out_act-self.rs_min[i]
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

    # Classify rings into categories
    def ring_classify(self):
        is_f=np.where(self.rs_min<=self.r_ds)[0]
        is_h=np.where(self.rs_min<=self.r_sh)[0]
        is_s=np.where((self.r_sh<self.rs_min) & (self.rs_min<=self.r_ds))[0]
        is_d=np.where(self.r_ds<self.rs_min)[0]
        # Central radius
        self.rs_f=self.rs[is_f]
        self.rs_h=self.rs[is_h]
        self.rs_s=self.rs[is_s]
        self.rs_d=self.rs[is_d]
        # Width (necessary to calculate weight)
        self.wids_f=self.wids[is_f]
        self.wids_h=self.wids[is_h]
        self.wids_s=self.wids[is_s]
        self.wids_d=self.wids[is_d]

    # Generator time-scale / Propagation time-scale
    def variability_frequency_calc(self, cb_f, m_f, cbp_f, mp_f, cb_d, m_d, cbp_d, mp_d):
        # Generator time-scale
        fs_vis_f=f_vis_calc(r=self.rs_f, lb=cb_f, m=m_f) #[c/Rg]
        fs_vis_d=f_vis_calc(r=self.rs_d, lb=cb_d, m=m_d) #[c/Rg]
        self.fs_vis=np.append(fs_vis_f, fs_vis_d)
        # Propagation time-scale
        fs_prop_f=f_vis_calc(r=self.rs_f, lb=cbp_f, m=mp_f) #[c/Rg]
        fs_prop_d=f_vis_calc(r=self.rs_d, lb=cbp_d, m=mp_d) #[c/Rg]
        self.fs_prop=np.append(fs_prop_f, fs_prop_d)
        # Propagation speed
        self.vs_prop=self.rs*self.fs_prop #[c]

    # Damping for all frequencies
    # Unnecessary! Will be removed in a future update.
    def damping_calc(self):
        self.cds=np.ones(self.n_ring)

    # Time taken for spectra to respond to mass accretion rate fluctuations
    # Unnecessary! Will be removed in a future update.
    def spec_lag_calc(self):
        self.lags=0.*np.ones(self.n_ring) #[s]
        self.lagrs=0.*np.ones(self.n_ring) #[s]
        self.lags_r=0.*np.ones(self.n_ring) #[s]
        self.lagrs_r=0.*np.ones(self.n_ring) #[s]

    # Intrinsic variability
    def variability_amplitude_calc(self, cf_var_f, dr_var_f, cf_var_d, dr_var_d):
        cfs_var_f=gauss(x=self.rs_f, norm=cf_var_f, mu=self.r_ds, sigma=dr_var_f)
        cfs_var_d=gauss(x=self.rs_d, norm=cf_var_d, mu=self.r_ds, sigma=dr_var_d)
        self.lfs_var=np.append(cfs_var_f, cfs_var_d)

    # Intrinsic variability
    def variability_amplitude_set(self, cf_var_d, cf_var_s, cf_var_h):
        cfs_var_h=cf_var_h*np.ones(len(self.rs_h))
        cfs_var_s=cf_var_s*np.ones(len(self.rs_s))
        cfs_var_d=cf_var_d*np.ones(len(self.rs_d))
        cfs_var=cfs_var_h
        cfs_var=np.append(cfs_var, cfs_var_s)
        cfs_var=np.append(cfs_var, cfs_var_d)
        self.lfs_var=cfs_var

    # Parameters of the top-hat impulse response for reverberation
    def imp_resp_set(self):
        # Delay
        t0s_d=self.t0_d*np.ones(len(self.rs_d))
        t0s_s=self.t0_s*np.ones(len(self.rs_s))
        t0s_h=self.t0_h*np.ones(len(self.rs_h))
        t0s=t0s_h
        t0s=np.append(t0s, t0s_s)
        t0s=np.append(t0s, t0s_d)
        self.t0s=t0s

        # Duration
        dt0s_d=self.dt0_d*np.ones(len(self.rs_d))
        dt0s_s=self.dt0_s*np.ones(len(self.rs_s))
        dt0s_h=self.dt0_h*np.ones(len(self.rs_h))
        dt0s=dt0s_h
        dt0s=np.append(dt0s, dt0s_s)
        dt0s=np.append(dt0s, dt0s_d)
        self.dt0s=dt0s

    # Final operation!
    def out2in(self):
        self.rs     =self.rs[::-1]
        self.rs_min =self.rs_min[::-1]
        self.wids   =self.wids[::-1]
        self.drs    =self.drs[::-1]
        self.rs_d   =self.rs_d[::-1]
        self.rs_f   =self.rs_f[::-1]
        self.rs_s   =self.rs_s[::-1]
        self.rs_h   =self.rs_h[::-1]
        self.fs_vis =self.fs_vis[::-1]
        self.fs_prop=self.fs_prop[::-1]
        self.vs_prop=self.vs_prop[::-1]
        self.eps    =self.eps[::-1]
        self.cds    =self.cds[::-1]
        self.lags   =self.lags[::-1]
        self.lagrs  =self.lagrs[::-1]
        self.lags_r =self.lags_r[::-1]
        self.lagrs_r=self.lagrs_r[::-1]
        self.lfs_var=self.lfs_var[::-1]
        self.t0s    =self.t0s[::-1]
        self.dt0s   =self.dt0s[::-1]
        self.ws     =self.ws[::-1]  # Subject band, Direct component
        #self.wrs    =self.wrs[::-1] # Subject band, Reflected component
        self.wrss   =self.wrss[::-1] # Subject band, Reflected component
        self.ws_r   =self.ws_r[::-1]  # Reference band, Direct component
        #self.wrs_r  =self.wrs_r[::-1] # Reference band, Reflected component
        self.wrss_r =self.wrss_r[::-1] # Reference band, Reflected component

# ---------------------------- #
# ----- Calculate weight ----- #
# ---------------------------- #
def preflow_weight_calc(\
    cs_d, cs_s, cs_h, cs_sr, cs_hr,\
    rs_d, rs_s, rs_h, wids_d, wids_s, wids_h,\
    stress, r_min, index_d, index_f):

    # --- Sensitivity x Spectrum --- #
    flux_d =cs_d
    flux_s =cs_s
    flux_h =cs_h
    flux_dr=0.
    flux_sr=cs_sr
    flux_hr=cs_hr

    # --- Calculate weight --- #
    # Stressed
    if stress==1:
        stress=True 
    # Stress-free
    elif stress==2:
        stress=False

    # Variable disk
    ws_d, ws_dr=weight_calc_final(\
        rs=rs_d, wids=wids_d, stress=stress, gamma=index_d, r_min=r_min,\
        flux=flux_d, flux_r=flux_dr)
    # Soft Compton
    ws_s, ws_sr=weight_calc_final(\
        rs=rs_s, wids=wids_s, stress=stress, gamma=index_f, r_min=r_min,\
        flux=flux_s, flux_r=flux_sr)
    # Hard Compton
    ws_h, ws_hr=weight_calc_final(\
        rs=rs_h, wids=wids_h, stress=stress, gamma=index_f, r_min=r_min,\
        flux=flux_h, flux_r=flux_hr)

    # Direct
    ws=ws_h
    ws=np.append(ws, ws_s)
    ws=np.append(ws, ws_d)
    # Reprocessed
    wrs=ws_hr
    wrs=np.append(wrs, ws_sr)
    wrs=np.append(wrs, ws_dr)

    return ws, wrs

###############
### PREFLOW ###
###############
def preflow_weight_calc(\
    cs_d, cs_s, cs_h,\
    csp_s, csp_h,\
    rs_d, rs_s, rs_h,\
    wids_d, wids_s, wids_h,\
    stress, r_min, index_d, index_f):

    # --- Calculate weight --- #
    # Stressed
    if stress==1:
        stress=True 
    # Stress-free
    elif stress==2:
        stress=False

    # Variable disk
    ws_d=weight_dir_calc_final(\
        rs=rs_d, wids=wids_d, stress=stress, gamma=index_d, r_min=r_min,\
        flux=cs_d)
    # Soft Compton
    ws_s=weight_dir_calc_final(\
        rs=rs_s, wids=wids_s, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_s)
    # Hard Compton
    ws_h=weight_dir_calc_final(\
        rs=rs_h, wids=wids_h, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_h)

    # --- The case the seed photon variability does NOT affect direct emission ------- #
    if [csp_s, csp_h]==[0., 0.] or \
       len(rs_d)==0:
        pass
    # --- The case the seed photon variability DOES affect direct emission --- #
    else:
        csp_d =0.
        # --- Incorporate seed photon variability into the weight --- #
        ws_d[0]+=(csp_d+csp_s+csp_h) 

    # Direct
    ws=ws_h
    ws=np.append(ws, ws_s)
    ws=np.append(ws, ws_d)

    return ws

def preflow_ref_weight_calc(\
    cs_d, cs_s, cs_h, cs_sr, cs_hr,\
    csp_s, csp_h, csp_sr, csp_hr,\
    rs_d,    rs_s,    rs_h,      wids_d,  wids_s,   wids_h,\
    stress,  r_min,   index_d,   index_f,\
    t0_d,    dt0_d,   t0_s,      dt0_s,   t0_h,     dt0_h,\
    fs):

    # --- Calculate weight --- #
    cs_dr=0.
    # Stressed
    if stress==1:
        stress=True 
    # Stress-free
    elif stress==2:
        stress=False

    # Variable disk
    ws_d, wss_dr=weight_dir_ref_calc_final(\
        rs=rs_d, wids=wids_d, stress=stress, gamma=index_d, r_min=r_min,\
        flux=cs_d, flux_r=cs_dr,\
        t0=t0_d, dt0=dt0_d,\
        fs=fs)
    # Soft Compton
    ws_s, wss_sr=weight_dir_ref_calc_final(\
        rs=rs_s, wids=wids_s, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_s, flux_r=cs_sr,\
        t0=t0_s, dt0=dt0_s,\
        fs=fs)
    # Hard Compton
    ws_h, wss_hr=weight_dir_ref_calc_final(\
        rs=rs_h, wids=wids_h, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_h, flux_r=cs_hr,\
        t0=t0_h, dt0=dt0_h,\
        fs=fs)

    # --- The case the seed photon variability does NOT affect direct emission ------- #
    if [csp_s, csp_h]==[0., 0., 0., 0.] or \
       len(rs_d)==0:
        pass
    # --- The case the seed photon variability DOES affect direct emission --- #
    else:
        csp_d =0.
        # --- Incorporate seed photon variability into the weight --- #
        ws_d[0]+=(csp_d+csp_s+csp_h) 

    # --- The case the seed photon variability does NOT affect reverberation ------- #
    #     1. Reverberation does not exist.                                            #
    #     2. Reverberation is completely insensitive to the seed photon variability.  #
    #     3. Variable disk does not exist.                                            #
    # ------------------------------------------------------------------------------- #
    if [csp_sr, csp_hr]==[0., 0.] or \
       len(rs_d)==0:
           pass
    # --- The case the seed photon variability DOES affect reverberation --- #
    else:
        # --- Incorporate seed photon variability into the weight --- #
        # Frequency-dependent weight for the reflection related to the seed photon variability
        ws_tilde_rs=csp_sr*ch_rep_calc(f=fs, t0=t0_s, dt0=dt0_s)+csp_hr*ch_rep_calc(f=fs, t0=t0_h, dt0=dt0_h)
        wss_dr[0]+=ws_tilde_rs

    # Direct
    ws=ws_h
    ws=np.append(ws, ws_s)
    ws=np.append(ws, ws_d)
    # Reprocessed
    wrss=wss_hr
    wrss=np.vstack((wrss, wss_sr))
    wrss=np.vstack((wrss, wss_dr))

    #print(ws_h)
    #print(ws_s)
    #print(ws_d)

    return ws, wrss

################
### PREFLOWS ###
################
def preflows_weight_calc(\
    e_min, e_max, pars_spec,\
    eta0_d,  eta1_d,  eta0_s,  eta1_s, eta0_h, eta1_h,\
    etap0_s, etap1_s, etap0_h, etap1_h,\
    rs_d, rs_s, rs_h, wids_d,  wids_s, wids_h,\
    stress, r_min, index_d, index_f):

    specx=SpectralModelXspec()

    # --- Calculate flux --- #
    es=np.array([e_min, e_max])
    de=e_max-e_min

    ktbbd =pars_spec[0] 
    normd =pars_spec[1] 
    gammas=pars_spec[2] 
    ecuts =pars_spec[3] 
    norms =pars_spec[4] 
    gammah=pars_spec[5] 
    ecuth =pars_spec[6] 
    normh =pars_spec[7] 
    incl  =pars_spec[8]
    a     =pars_spec[9]
    afe   =pars_spec[10]
    rins  =pars_spec[11]
    routs =pars_spec[12]
    indexs=pars_spec[13]
    logxis=pars_spec[14]
    normsr=pars_spec[15]
    rinh  =pars_spec[16]
    routh =pars_spec[17]
    indexh=pars_spec[18]
    logxih=pars_spec[19]
    normhr=pars_spec[20]

    flux_d =specx.diskbb_spec   (es=es, temp=ktbbd, norm=normd)[0]/de
    flux_s =specx.cutoffpl_spec (es=es, gamma=gammas, ecut=ecuts, norm=norms)[0]/de
    flux_h =specx.cutoffpl_spec (es=es, gamma=gammah, ecut=ecuth, norm=normh)[0]/de

    # Normalize such that total corresponds to unity.
    flux_tot=flux_d+flux_s+flux_h
    flux_d /=flux_tot
    flux_s /=flux_tot
    flux_h /=flux_tot

    # --- Calculate sensitivity parameter --- #
    # Averaged \eta (sensitivity parameter)
    # No upper and lower bounds
    eta_d =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_d, eta1=eta1_d)
    eta_s =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_s, eta1=eta1_s)
    eta_h =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_h, eta1=eta1_h)

    # --- Sensitivity x Spectrum --- #
    cs_d =flux_d *eta_d
    cs_s =flux_s *eta_s
    cs_h =flux_h *eta_h

    # --- Calculate weight --- #
    # Stressed
    if stress==1:
        stress=True 
    # Stress-free
    elif stress==2:
        stress=False

    # Variable disk
    ws_d=weight_dir_calc_final(\
        rs=rs_d, wids=wids_d, stress=stress, gamma=index_d, r_min=r_min,\
        flux=cs_d)
    # Soft Compton
    ws_s=weight_dir_calc_final(\
        rs=rs_s, wids=wids_s, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_s)
    # Hard Compton
    ws_h=weight_dir_calc_final(\
        rs=rs_h, wids=wids_h, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_h)

    # --- The case the seed photon variability does NOT affect direct emission ------- #
    if [etap0_s, etap1_s, etap0_h, etap1_h]==[0., 0., 0., 0.] or \
       len(rs_d)==0:
        pass
    # --- The case the seed photon variability DOES affect direct emission --- #
    else:
        etap_d =0. # Variable disk does not react to the seed photon variability
        etap_s =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_s, eta1=etap1_s)
        etap_h =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_h, eta1=etap1_h)
    
        csp_d =flux_d *etap_d # 0.
        csp_s =flux_s *etap_s
        csp_h =flux_h *etap_h

        # --- Incorporate seed photon variability into the weight --- #
        ws_d[0]+=(csp_d+csp_s+csp_h) 

    #print(ws_h)
    #print(ws_s)
    #print(ws_d)

    # Direct
    ws=ws_h
    ws=np.append(ws, ws_s)
    ws=np.append(ws, ws_d)

    return ws

def preflows_ref_weight_calc(\
    e_min,   e_max,   pars_spec,\
    eta0_d,  eta1_d,  eta0_s,    eta1_s,  eta0_h,   eta1_h,   eta0_sr,  eta1_sr, eta0_hr, eta1_hr,\
    etap0_s, etap1_s, etap0_h,   etap1_h, etap0_sr, etap1_sr, etap0_hr, etap1_hr,\
    rs_d,    rs_s,    rs_h,      wids_d,  wids_s,   wids_h,\
    stress,  r_min,   index_d,   index_f,\
    t0_d,    dt0_d,   t0_s,      dt0_s,   t0_h,     dt0_h,\
    fs):

    specx=SpectralModelXspec()

    # --- Calculate flux --- #
    es=np.array([e_min, e_max])
    de=e_max-e_min

    ktbbd =pars_spec[0] 
    normd =pars_spec[1] 
    gammas=pars_spec[2] 
    ecuts =pars_spec[3] 
    norms =pars_spec[4] 
    gammah=pars_spec[5] 
    ecuth =pars_spec[6] 
    normh =pars_spec[7] 
    incl  =pars_spec[8]
    a     =pars_spec[9]
    afe   =pars_spec[10]
    rins  =pars_spec[11]
    routs =pars_spec[12]
    indexs=pars_spec[13]
    logxis=pars_spec[14]
    normsr=pars_spec[15]
    rinh  =pars_spec[16]
    routh =pars_spec[17]
    indexh=pars_spec[18]
    logxih=pars_spec[19]
    normhr=pars_spec[20]

    flux_d =specx.diskbb_spec   (es=es, temp=ktbbd, norm=normd)[0]/de
    flux_s =specx.cutoffpl_spec (es=es, gamma=gammas, ecut=ecuts, norm=norms)[0]/de
    flux_h =specx.cutoffpl_spec (es=es, gamma=gammah, ecut=ecuth, norm=normh)[0]/de
    flux_dr=0.
    flux_sr=specx.relxill_spec(\
        es=es,        incl=incl,    a=a,\
        r_in=rins,    r_out=routs,  index=indexs,\
        gamma=gammas, logxi=logxis,\
        cafe=afe,     kte=ecuts,    norm=normsr)[0]/de
    flux_hr=specx.relxill_spec(\
        es=es,        incl=incl,    a=a,\
        r_in=rinh,    r_out=routh,  index=indexh,\
        gamma=gammah, logxi=logxih,\
        cafe=afe,     kte=ecuth,    norm=normhr)[0]/de

    # Normalize such that total corresponds to unity.
    flux_tot=flux_d+flux_s+flux_h+flux_dr+flux_sr+flux_hr
    flux_d /=flux_tot
    flux_s /=flux_tot
    flux_h /=flux_tot
    flux_dr/=flux_tot
    flux_sr/=flux_tot
    flux_hr/=flux_tot

    # --- Calculate sensitivity parameter --- #
    ### Mdot fluctuations ###
    # Averaged \eta (sensitivity parameter)
    # No upper and lower bounds
    eta_d =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_d, eta1=eta1_d)
    eta_s =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_s, eta1=eta1_s)
    eta_h =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_h, eta1=eta1_h)
    eta_dr=1.
    eta_sr=calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_sr, eta1=eta1_sr)
    eta_hr=calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_hr, eta1=eta1_hr)

    # --- Sensitivity x Spectrum --- #
    ### Mdot fluctuations ###
    cs_d =flux_d *eta_d
    cs_s =flux_s *eta_s
    cs_h =flux_h *eta_h
    cs_dr=flux_dr*eta_dr # 0.
    cs_sr=flux_sr*eta_sr
    cs_hr=flux_hr*eta_hr

    # --- Calculate weight --- #
    # Stressed
    if stress==1:
        stress=True 
    # Stress-free
    elif stress==2:
        stress=False

    # Variable disk
    ws_d, wss_dr=weight_dir_ref_calc_final(\
        rs=rs_d, wids=wids_d, stress=stress, gamma=index_d, r_min=r_min,\
        flux=cs_d, flux_r=cs_dr,\
        t0=t0_d, dt0=dt0_d,\
        fs=fs)
    # Soft Compton
    ws_s, wss_sr=weight_dir_ref_calc_final(\
        rs=rs_s, wids=wids_s, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_s, flux_r=cs_sr,\
        t0=t0_s, dt0=dt0_s,\
        fs=fs)
    # Hard Compton
    ws_h, wss_hr=weight_dir_ref_calc_final(\
        rs=rs_h, wids=wids_h, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_h, flux_r=cs_hr,\
        t0=t0_h, dt0=dt0_h,\
        fs=fs)

    # --- The case the seed photon variability does NOT affect direct emission ------- #
    if [etap0_s, etap1_s, etap0_h, etap1_h]==[0., 0., 0., 0.] or \
       len(rs_d)==0:
        pass
    # --- The case the seed photon variability DOES affect direct emission --- #
    else:
        etap_d =0. # Variable disk does not react to the seed photon variability
        etap_s =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_s, eta1=etap1_s)
        etap_h =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_h, eta1=etap1_h)
    
        csp_d =flux_d *etap_d # 0.
        csp_s =flux_s *etap_s
        csp_h =flux_h *etap_h

        # --- Incorporate seed photon variability into the weight --- #
        ws_d[0]+=(csp_d+csp_s+csp_h) 

    # --- The case the seed photon variability does NOT affect reverberation ------- #
    #     1. Reverberation does not exist.                                            #
    #     2. Reverberation is completely insensitive to the seed photon variability.  #
    #     3. Variable disk does not exist.                                            #
    # ------------------------------------------------------------------------------- #
    if [normsr, normhr]==[0., 0.] or \
       [etap0_sr, etap1_sr, etap0_hr, etap1_hr]==[0., 0., 0., 0.] or \
       len(rs_d)==0:
           pass
    # --- The case the seed photon variability DOES affect reverberation --- #
    else:
        ### Seed photon fluctuations ###
        etap_dr=1.
        etap_sr=calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_sr, eta1=etap1_sr)
        etap_hr=calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_hr, eta1=etap1_hr)

        csp_dr=flux_dr*etap_dr # 0.
        csp_sr=flux_sr*etap_sr
        csp_hr=flux_hr*etap_hr

        # --- Incorporate seed photon variability into the weight --- #
        # Frequency-dependent weight for the reflection related to the seed photon variability
        ws_tilde_rs=csp_sr*ch_rep_calc(f=fs, t0=t0_s, dt0=dt0_s)+csp_hr*ch_rep_calc(f=fs, t0=t0_h, dt0=dt0_h)
        wss_dr[0]+=ws_tilde_rs

    # Direct
    ws=ws_h
    ws=np.append(ws, ws_s)
    ws=np.append(ws, ws_d)
    # Reprocessed
    wrss=wss_hr
    wrss=np.vstack((wrss, wss_sr))
    wrss=np.vstack((wrss, wss_dr))

    #print(ws_h)
    #print(ws_s)
    #print(ws_d)

    return ws, wrss

##################
### PREFLOWSCP ###
##################
def preflowscp_weight_calc(\
    e_min, e_max, pars_spec,\
    eta0_d,  eta1_d,  eta0_s,  eta1_s, eta0_h, eta1_h,\
    etap0_s, etap1_s, etap0_h, etap1_h,\
    rs_d, rs_s, rs_h, wids_d,  wids_s, wids_h,\
    stress, r_min, index_d, index_f):

    specx=SpectralModelXspec()

    # --- Calculate flux --- #
    es=np.array([e_min, e_max])
    de=e_max-e_min

    ktbbd =pars_spec[0] 
    normd =pars_spec[1] 
    gammas=pars_spec[2] 
    ktbbs =pars_spec[3] 
    ktes  =pars_spec[4] 
    norms =pars_spec[5] 
    gammah=pars_spec[6] 
    ktbbh =pars_spec[7] 
    kteh  =pars_spec[8] 
    normh =pars_spec[9] 
    incl  =pars_spec[10]
    a     =pars_spec[11]
    afe   =pars_spec[12]
    rins  =pars_spec[13]
    routs =pars_spec[14]
    indexs=pars_spec[15]
    logxis=pars_spec[16]
    logns =pars_spec[17]
    normsr=pars_spec[18]
    rinh  =pars_spec[19]
    routh =pars_spec[20]
    indexh=pars_spec[21]
    logxih=pars_spec[22]
    lognh =pars_spec[23]
    normhr=pars_spec[24]

    flux_d =specx.diskbb_spec   (es=es, temp=ktbbd, norm=normd)[0]/de
    flux_s =specx.nthcomp_spec  (\
        es=es, gamma=gammas, ktbb=ktbbs, kte=ktes, norm=norms)[0]/de
    flux_h =specx.nthcomp_spec  (\
        es=es, gamma=gammah, ktbb=ktbbh, kte=ktes, norm=normh)[0]/de

    # Normalize such that total corresponds to unity.
    flux_tot=flux_d+flux_s+flux_h
    flux_d /=flux_tot
    flux_s /=flux_tot
    flux_h /=flux_tot

    # --- Calculate sensitivity parameter --- #
    # Averaged \eta (sensitivity parameter)
    # No upper and lower bounds
    eta_d =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_d, eta1=eta1_d)
    eta_s =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_s, eta1=eta1_s)
    eta_h =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_h, eta1=eta1_h)

    # --- Sensitivity x Spectrum --- #
    cs_d =flux_d *eta_d
    cs_s =flux_s *eta_s
    cs_h =flux_h *eta_h

    # --- Calculate weight --- #
    # Stressed
    if stress==1:
        stress=True 
    # Stress-free
    elif stress==2:
        stress=False

    # Variable disk
    ws_d=weight_dir_calc_final(\
        rs=rs_d, wids=wids_d, stress=stress, gamma=index_d, r_min=r_min,\
        flux=cs_d)
    # Soft Compton
    ws_s=weight_dir_calc_final(\
        rs=rs_s, wids=wids_s, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_s)
    # Hard Compton
    ws_h=weight_dir_calc_final(\
        rs=rs_h, wids=wids_h, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_h)

    # --- The case the seed photon variability does NOT affect direct emission ------- #
    if [etap0_s, etap1_s, etap0_h, etap1_h]==[0., 0., 0., 0.] or \
       len(rs_d)==0:
        pass
    # --- The case the seed photon variability DOES affect direct emission --- #
    else:
        etap_d =0. # Variable disk does not react to the seed photon variability
        etap_s =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_s, eta1=etap1_s)
        etap_h =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_h, eta1=etap1_h)
    
        csp_d =flux_d *etap_d # 0.
        csp_s =flux_s *etap_s
        csp_h =flux_h *etap_h

        # --- Incorporate seed photon variability into the weight --- #
        ws_d[0]+=(csp_d+csp_s+csp_h) 

    #print(ws_h)
    #print(ws_s)
    #print(ws_d)

    # Direct
    ws=ws_h
    ws=np.append(ws, ws_s)
    ws=np.append(ws, ws_d)

    return ws

def preflowscp_ref_weight_calc(\
    e_min,   e_max,   pars_spec,\
    eta0_d,  eta1_d,  eta0_s,    eta1_s,  eta0_h,   eta1_h,   eta0_sr,  eta1_sr, eta0_hr, eta1_hr,\
    etap0_s, etap1_s, etap0_h,   etap1_h, etap0_sr, etap1_sr, etap0_hr, etap1_hr,\
    rs_d,    rs_s,    rs_h,      wids_d,  wids_s,   wids_h,\
    stress,  r_min,   index_d,   index_f,\
    t0_d,    dt0_d,   t0_s,      dt0_s,   t0_h,     dt0_h,\
    fs):

    specx=SpectralModelXspec()

    # --- Calculate flux --- #
    es=np.array([e_min, e_max])
    de=e_max-e_min

    ktbbd =pars_spec[0] 
    normd =pars_spec[1] 
    gammas=pars_spec[2] 
    ktbbs =pars_spec[3] 
    ktes  =pars_spec[4] 
    norms =pars_spec[5] 
    gammah=pars_spec[6] 
    ktbbh =pars_spec[7] 
    kteh  =pars_spec[8] 
    normh =pars_spec[9] 
    incl  =pars_spec[10]
    a     =pars_spec[11]
    afe   =pars_spec[12]
    rins  =pars_spec[13]
    routs =pars_spec[14]
    indexs=pars_spec[15]
    logxis=pars_spec[16]
    logns =pars_spec[17]
    normsr=pars_spec[18]
    rinh  =pars_spec[19]
    routh =pars_spec[20]
    indexh=pars_spec[21]
    logxih=pars_spec[22]
    lognh =pars_spec[23]
    normhr=pars_spec[24]

    flux_d =specx.diskbb_spec   (es=es, temp=ktbbd, norm=normd)[0]/de
    flux_s =specx.nthcomp_spec  (\
        es=es, gamma=gammas, ktbb=ktbbs, kte=ktes, norm=norms)[0]/de
    flux_h =specx.nthcomp_spec  (\
        es=es, gamma=gammah, ktbb=ktbbh, kte=ktes, norm=normh)[0]/de
    flux_dr=0. # No reflection associated with the variable disk
    flux_sr=specx.relxillcp_spec(\
        es=es,        incl=incl,    a=a,\
        r_in=rins,    r_out=routs,  index=indexs,\
        gamma=gammas, logxi=logxis, logcn=logns,\
        cafe=afe,     kte=ktes,     norm=normsr)[0]/de
    flux_hr=specx.relxillcp_spec(\
        es=es,        incl=incl,    a=a,\
        r_in=rinh,    r_out=routh,  index=indexh,\
        gamma=gammah, logxi=logxih, logcn=lognh,\
        cafe=afe,     kte=kteh,     norm=normhr)[0]/de

    # Normalize such that total corresponds to unity.
    flux_tot=flux_d+flux_s+flux_h+flux_dr+flux_sr+flux_hr
    flux_d /=flux_tot
    flux_s /=flux_tot
    flux_h /=flux_tot
    flux_dr/=flux_tot
    flux_sr/=flux_tot
    flux_hr/=flux_tot

    # --- Calculate sensitivity parameter --- #
    ### Mdot fluctuations ###
    # Averaged \eta (sensitivity parameter)
    # No upper and lower bounds
    eta_d =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_d, eta1=eta1_d)
    eta_s =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_s, eta1=eta1_s)
    eta_h =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_h, eta1=eta1_h)
    eta_dr=1.
    eta_sr=calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_sr, eta1=eta1_sr)
    eta_hr=calc_eta_ave(e_min=e_min, e_max=e_max, eta0=eta0_hr, eta1=eta1_hr)

    # --- Sensitivity x Spectrum --- #
    ### Mdot fluctuations ###
    cs_d =flux_d *eta_d
    cs_s =flux_s *eta_s
    cs_h =flux_h *eta_h
    cs_dr=flux_dr*eta_dr # 0.
    cs_sr=flux_sr*eta_sr
    cs_hr=flux_hr*eta_hr

    # --- Calculate weight --- #
    # Stressed
    if stress==1:
        stress=True 
    # Stress-free
    elif stress==2:
        stress=False

    # Variable disk
    ws_d, wss_dr=weight_dir_ref_calc_final(\
        rs=rs_d, wids=wids_d, stress=stress, gamma=index_d, r_min=r_min,\
        flux=cs_d, flux_r=cs_dr,\
        t0=t0_d, dt0=dt0_d,\
        fs=fs)
    # Soft Compton
    ws_s, wss_sr=weight_dir_ref_calc_final(\
        rs=rs_s, wids=wids_s, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_s, flux_r=cs_sr,\
        t0=t0_s, dt0=dt0_s,\
        fs=fs)
    # Hard Compton
    ws_h, wss_hr=weight_dir_ref_calc_final(\
        rs=rs_h, wids=wids_h, stress=stress, gamma=index_f, r_min=r_min,\
        flux=cs_h, flux_r=cs_hr,\
        t0=t0_h, dt0=dt0_h,\
        fs=fs)

    # --- The case the seed photon variability does NOT affect direct emission ------- #
    if [etap0_s, etap1_s, etap0_h, etap1_h]==[0., 0., 0., 0.] or \
       len(rs_d)==0:
        pass
    # --- The case the seed photon variability DOES affect direct emission --- #
    else:
        etap_d =0. # Variable disk does not react to the seed photon variability
        etap_s =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_s, eta1=etap1_s)
        etap_h =calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_h, eta1=etap1_h)
    
        csp_d =flux_d *etap_d # 0.
        csp_s =flux_s *etap_s
        csp_h =flux_h *etap_h

        # --- Incorporate seed photon variability into the weight --- #
        ws_d[0]+=(csp_d+csp_s+csp_h) 

    # --- The case the seed photon variability does NOT affect reverberation ------- #
    #     1. Reverberation does not exist.                                            #
    #     2. Reverberation is completely insensitive to the seed photon variability.  #
    #     3. Variable disk does not exist.                                            #
    # ------------------------------------------------------------------------------- #
    if [normsr, normhr]==[0., 0.] or \
       [etap0_sr, etap1_sr, etap0_hr, etap1_hr]==[0., 0., 0., 0.] or \
       len(rs_d)==0:
           pass
    # --- The case the seed photon variability DOES affect reverberation --- #
    else:
        ### Seed photon fluctuations ###
        etap_dr=1.
        etap_sr=calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_sr, eta1=etap1_sr)
        etap_hr=calc_eta_ave(e_min=e_min, e_max=e_max, eta0=etap0_hr, eta1=etap1_hr)

        csp_dr=flux_dr*etap_dr # 0.
        csp_sr=flux_sr*etap_sr
        csp_hr=flux_hr*etap_hr

        # --- Incorporate seed photon variability into the weight --- #
        # Frequency-dependent weight for the reflection related to the seed photon variability
        ws_tilde_rs=csp_sr*ch_rep_calc(f=fs, t0=t0_s, dt0=dt0_s)+csp_hr*ch_rep_calc(f=fs, t0=t0_h, dt0=dt0_h)
        wss_dr[0]+=ws_tilde_rs

    # Direct
    ws=ws_h
    ws=np.append(ws, ws_s)
    ws=np.append(ws, ws_d)
    # Reprocessed
    wrss=wss_hr
    wrss=np.vstack((wrss, wss_sr))
    wrss=np.vstack((wrss, wss_dr))

    #print(ws_h)
    #print(ws_s)
    #print(ws_d)

    return ws, wrss

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

def gauss(x, norm, mu, sigma):
    y=norm*np.exp(-((x-mu)**2)/(2.*(sigma**2)))
    return y

def weight_calc_final(rs, wids, stress, gamma, r_min, flux, flux_r):
    es_dis=epsilon_calc(r=rs, stress=stress, gamma=gamma, r_min=r_min)*2.*np.pi*rs*wids
    es_dis_tot=np.sum(es_dis) # Total energy dissipated (normalization)
    # Direct
    ws  =flux  *es_dis/es_dis_tot
    # Reprocessed
    ws_r=flux_r*es_dis/es_dis_tot
    return ws, ws_r

def weight_dir_calc_final(rs, wids, stress, gamma, r_min, flux):
    es_dis=epsilon_calc(r=rs, stress=stress, gamma=gamma, r_min=r_min)*2.*np.pi*rs*wids
    es_dis_tot=np.sum(es_dis) # Total energy dissipated (normalization)
    # Direct
    ws  =flux  *es_dis/es_dis_tot
    return ws

def weight_dir_ref_calc_final(\
    rs, wids, stress, gamma, r_min,\
    flux, flux_r,\
    t0, dt0,\
    fs):
    es_dis=epsilon_calc(r=rs, stress=stress, gamma=gamma, r_min=r_min)*2.*np.pi*rs*wids
    es_dis_tot=np.sum(es_dis) # Total energy dissipated (normalization)
    # Direct
    ws  =flux  *es_dis/es_dis_tot
    # Reprocessed
    ws_r=flux_r*es_dis/es_dis_tot
    first=True
    for w_r in ws_r:
        ws_r_f=w_r*ch_rep_calc(f=fs, t0=t0, dt0=dt0)
        if first==True:
            first=False
            wss_r=ws_r_f
        else:
            wss_r=np.vstack((wss_r, ws_r_f))
    return ws, wss_r

# ------------------------------------ #
# ----- Mass accretion rate PSD ------ #
# ------------------------------------ #
class FluPro:
    def __init__(self):
        self.mu=1. # fixed without the loss of generarity
        self.sigma=0
        self.norm_psd=0.
        self.n_data=0
        self.dt=0
        self.f_set_done=False
        self.f_vis_set_done=False
        self.sigma_set_done=False
        self.psd_wo_prop_done=False

    def f_set(self, fs):
        self.fs=fs #[Hz]
        self.f_set_done=True

    def f_vis_set(self, fs_vis):
        self.fs_vis=fs_vis #[c_rg]
        self.f_vis_set_done=True

    def sigma_set(self, sigs):
        self.sigs=sigs 
        self.sigma_set_done=True

    # Calculate PSD without propagation
    def psd_wo_prop(self, c_rg):
        if self.f_set_done==False:
            print('Error: parameter is not set.')
            sys.exit()
        if self.f_vis_set_done==False:
            print('Error: parameter is not set.')
            sys.exit()
        if self.sigma_set_done==False:
            print('Error: parameter is not set.')
            sys.exit()
        for i, (f_vis, sig) in enumerate(zip(self.fs_vis, self.sigs)):
            df=f_vis*c_rg # [Hz]
            # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2
            psd_intr=lorentz(f=self.fs, mu=self.mu, sigma=sig, f_c=0, df=df) 

            if i==0:
                self.psds_intr=psd_intr
            else:
                self.psds_intr=np.vstack((self.psds_intr, psd_intr))
        self.psd_wo_prop_done=True

    # Calculate PSD with propagation (No damping)
    #def psd_w_prop(self):
    #    if self.psd_wo_prop_done==False:
    #        print('Error: PSD without propagation is not calculated.')
    #        sys.exit()

    #    self.norm_psd=2.*self.dt/((self.mu**2)*self.n_data)
    #    for i in range(len(self.fs_vis)):
    #        b2s=(1./self.norm_psd)*self.psds_intr[i] # Modulus square of the Fourier transform 
    #        # ---------------------- #
    #        # --- Outermost ring --- #
    #        # ---------------------- #
    #        if i==0:
    #            #lm2s=self.mdot0*b2s
    #            lm2s=b2s # Modulus square of the Fourier transform

    #            # Without damping
    #            #psd_prop=self.norm_psd*lm2s # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2
    #            # With damping
    #            psd_prop=self.norm_psd*lm2s # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2

    #            self.psds_prop=psd_prop
    #            psd_prop_pre=psd_prop

    #        # ------------------- #
    #        # --- Inner rings --- #
    #        # ------------------- #
    #        else:
    #            c2s=(1./self.norm_psd)*psd_prop_pre # Modulus square of the Fourier transform
    #            fs_exte, b2s_exte=ft_mod2_unfold(fs=self.fs, b2s=b2s, mu=self.mu)
    #            fs_exte, c2s_exte=ft_mod2_unfold(fs=self.fs, b2s=c2s, mu=self.mu)

    #            lm2s_exte=conv_calc_eff(bs=b2s_exte, cs=c2s_exte)
    #            fs, lm2s=ft_mod2_fold(fs=fs_exte, bs=lm2s_exte)
    #            # Modulus square of the Fourier transform
    #            # |X_{j}|^2=( (|A_{j}|^2) \odot (|B_{j}|^2) ) / N^2 in our definition of Fourier transform
    #            lm2s/=self.n_data**2 

    #            # Without damping
    #            #psd_prop=self.norm_psd*lm2s # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2
    #            # With damping
    #            psd_prop=self.norm_psd*lm2s # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2

    #            self.psds_prop=np.vstack((self.psds_prop, psd_prop))
    #            psd_prop_pre=psd_prop

    ## !!! Wrong code !!!
    ### 2022/06/04 ###
    ### ------------------------------------------------------ ###
    ### Employ the green function of Rapisarda et al.2017a (3) ###
    ### ------------------------------------------------------ ###
    ## Calculate PSD with propagation
    #def psd_w_prop(self, cs, fs_prop, dr_r, rg_c):
    #    if self.psd_wo_prop_done==False:
    #        print('Error: PSD without propagation is not calculated.')
    #        sys.exit()

    #    self.norm_psd=2.*self.dt/((self.mu**2)*self.n_data)
    #    for i in range(len(self.fs_vis)):
    #        b2s=(1./self.norm_psd)*self.psds_intr[i] # Modulus square of the Fourier transform 
    #        # ---------------------- #
    #        # --- Outermost ring --- #
    #        # ---------------------- #
    #        if i==0:
    #            lm2s=b2s # Modulus square of the Fourier transform
    #            psd_prop=self.norm_psd*lm2s # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2
    #            self.psds_prop=psd_prop
    #            psd_prop_pre=psd_prop

    #        # ------------------- #
    #        # --- Inner rings --- #
    #        # ------------------- #
    #        else:
    #            c2s=(1./self.norm_psd)*psd_prop_pre # Modulus square of the Fourier transform
    #            fs_exte, b2s_exte=ft_mod2_unfold(fs=self.fs, b2s=b2s, mu=self.mu)
    #            fs_exte, c2s_exte=ft_mod2_unfold(fs=self.fs, b2s=c2s, mu=self.mu)

    #            lm2s_exte=conv_calc_eff(bs=b2s_exte, cs=c2s_exte)
    #            fs, lm2s=ft_mod2_fold(fs=fs_exte, bs=lm2s_exte)
    #            # Modulus square of the Fourier transform
    #            # |X_{j}|^2=( (|A_{j}|^2) \odot (|B_{j}|^2) ) / N^2 in our definition of Fourier transform
    #            lm2s/=self.n_data**2 

    #            # Propagation time from the previous ring to current ring = (dr/r)*(1/fvisc(current ring))
    #            t_prop=dr_r/fs_prop[i-1] #[Rg/c]
    #            t_prop*=rg_c #[s]
    #            # Green function |G(r_{n-1}, r_n, f)|^2
    #            cgs_mod2=green_function_ft_mod2(f=self.fs, cs=cs, t_prop=t_prop)

    #            psd_prop=self.norm_psd*lm2s*cgs_mod2 # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2

    #            self.psds_prop=np.vstack((self.psds_prop, psd_prop))
    #            psd_prop_pre=psd_prop

    ### 2022/08/23 ###
    ### ------------------------------------------------------ ###
    ### Employ the green function of Rapisarda et al.2017a (3) ###
    ### ------------------------------------------------------ ###
    # Calculate PSD with propagation
    def psd_w_prop(self, cs, fs_prop, dr_r, rg_c):
        if self.psd_wo_prop_done==False:
            print('Error: PSD without propagation is not calculated.')
            sys.exit()

        self.norm_psd=2.*self.dt/((self.mu**2)*self.n_data)
        for i in range(len(self.fs_vis)):
            b2s=(1./self.norm_psd)*self.psds_intr[i] # Modulus square of the Fourier transform 
            # ---------------------- #
            # --- Outermost ring --- #
            # ---------------------- #
            if i==0:
                lm2s=b2s # Modulus square of the Fourier transform
                psd_prop=self.norm_psd*lm2s # \int _{0} ^{\infty} df P(f)=(\sigma/\mu)^2
                self.psds_prop=psd_prop
                psd_prop_pre=psd_prop

            # ------------------- #
            # --- Inner rings --- #
            # ------------------- #
            else:
                c2s=(1./self.norm_psd)*psd_prop_pre # Modulus square of the Fourier transform
                fs_exte, b2s_exte=ft_mod2_unfold(fs=self.fs, b2s=b2s, mu=self.mu)
                #fs_exte, c2s_exte=ft_mod2_unfold(fs=self.fs, b2s=c2s, mu=self.mu)

                # Propagation time from the previous ring to current ring = (dr/r)*(1/fvisc(current ring))
                t_prop=dr_r/fs_prop[i-1] #[Rg/c]
                t_prop*=rg_c #[s]
                # Green function |G(r_{n-1}, r_n, f)|^2
                cgs_mod2=green_function_ft_mod2(f=self.fs, cs=cs, t_prop=t_prop)
                # Format for FFT
                fs_exte, c2s_exte=ft_mod2_unfold(fs=self.fs, b2s=c2s*cgs_mod2, mu=self.mu)

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

    #def psd_norm_set(self, dt, n_data):
    #    # PSD [(sigma/mean)^2/Hz]
    #    self.norm_psd=2.*dt/((self.mu_fl**2)*n_data) 

    def psd_norm_set(self):
        # PSD [(sigma/mean)^2/Hz]
        self.norm_psd=2.*self.dt/((self.mu_fl**2)*self.n_data) 

    #def csd_norm_set(self, dt, n_data):
    #    # CSD [(sigma/mean)^2/Hz]
    #    self.norm_csd=2.*dt/((self.mu_fl_ref*self.mu_fl)*n_data) 

    def csd_norm_set(self):
        # CSD [(sigma/mean)^2/Hz]
        self.norm_csd=2.*self.dt/((self.mu_fl_r*self.mu_fl)*self.n_data) 

    def norm_rep_set(self, dt0, w_flow_tot):
        #C(E) in the impulse response
        ref=self.speceff_sref+self.speceff_mref+self.speceff_href
        self.norm_rep=ref/(w_flow_tot*dt0) 

    def norm_rep_ref_set(self, dt0, w_flow_tot):
        #C(E) in the impulse response
        ref=self.speceff_sref_ref+self.speceff_mref_ref+self.speceff_href_ref
        self.norm_rep_ref=ref/(w_flow_tot*dt0) 

    # Calculation of PSD without reverberation
    def psd_flux_calc(self):
        if self.norm_psd==0:
            print('Error: normalization of PSD is not set.')
            sys.exit()
        lf2s=lf2_calc(\
            fs     =self.fs,\
            n_r    =self.n_ring,\
            ws     =self.ws,\
            lm2s   =self.lm2s,\
            fs_prop=self.fs_prop,\
            cs     =self.cs,\
            dr_r   =self.dr_r,\
            rg_c   =self.rg_c)
        # Normalize
        self.psd_fl=self.norm_psd*lf2s

    # Calculation of PSD with reverberation
    def psd_flux_ref_calc(self):
        if self.norm_psd==0:
            print('Error: normalization of PSD is not set.')
            sys.exit()
        lf2s=lf2_ref_calc(\
            fs     =self.fs,\
            n_r    =self.n_ring,\
            ws     =self.ws,\
            wrss   =self.wrss,\
            lm2s   =self.lm2s,\
            fs_prop=self.fs_prop,\
            cs     =self.cs,\
            dr_r   =self.dr_r,\
            rg_c   =self.rg_c)
        # Normalize
        self.psd_fl=self.norm_psd*lf2s

    # Calculation of CSD without reverberation
    def csd_flux_calc(self):
        if self.norm_csd==0:
            print('Error: normalization of CSD is not set.')
            sys.exit()
        lflfs=lflf_calc(\
            fs     =self.fs,\
            n_r    =self.n_ring,\
            ws     =self.ws,\
            ws_r   =self.ws_r,\
            lm2s   =self.lm2s,\
            fs_prop=self.fs_prop,\
            cs     =self.cs,\
            dr_r   =self.dr_r,\
            rg_c   =self.rg_c)
        # Normalize
        self.csd_fl=self.norm_csd*lflfs
        self.csd_flux_calc_done=True

    # Calculation of CSD with reverberation
    def csd_flux_ref_calc(self):
        if self.norm_csd==0:
            print('Error: normalization of CSD is not set.')
            sys.exit()
        lflfs=lflf_ref_calc(\
            fs     =self.fs,\
            n_r    =self.n_ring,\
            ws     =self.ws,\
            wrss   =self.wrss,\
            ws_r   =self.ws_r,\
            wrss_r =self.wrss_r,\
            lm2s   =self.lm2s,\
            fs_prop=self.fs_prop,\
            cs     =self.cs,\
            dr_r   =self.dr_r,\
            rg_c   =self.rg_c)
        # Normalize
        self.csd_fl=self.norm_csd*lflfs
        self.csd_flux_calc_done=True

    def csd_convert(self):
        # Re[CSD]
        if self.quant==2:
            self.csd_fl=np.real(self.csd_fl)
        # Im[CSD]
        elif self.quant==3:
            if self.invert==1:
                self.csd_fl=np.imag(self.csd_fl)
            elif self.invert==2:
                self.csd_fl=-np.imag(self.csd_fl)
        # |CSD|
        elif self.quant==4:
            self.csd_fl=np.abs(self.csd_fl)
        # Phase lag
        elif self.quant==5:
            self.lagf_calc()
            if self.invert==1:
                self.csd_fl=self.phi
            elif self.invert==2:
                self.csd_fl=-self.phi
        # Time lag
        elif self.quant==6:
            self.lagf_calc()
            if self.invert==1:
                self.csd_fl=self.tau
            elif self.invert==2:
                self.csd_fl=-self.tau

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

# --------------------------------------- #
# ---------- Sensitivity (eta) ---------- #
# --------------------------------------- #
def calc_eta(e, eta0, eta1):
    eta=eta0+eta1*np.log10(e)
    return eta

def calc_eta_ave(e_min, e_max, eta0, eta1):
    eta=(eta1/np.log(10.))*(np.log(e_max)-(e_min/e_max)*np.log(e_min))/(1.-(e_min/e_max))\
        + (eta0-eta1/np.log(10.))
    return eta

# ------------------------------------ #
# ---------- Green function ---------- #
# ------------------------------------ #
# G(rn, f): Same as (3) in Rapisarda et al. 2017a
def green_function_ft(f, cs, t_prop):
    cg=np.exp(-cs*f*t_prop)*np.exp(-1j*2.*np.pi*f*t_prop)
    return cg

# |G(rn, f)|: Practical use
def green_function_ft_amp(f, cs, t_prop):
    cg_amp=np.exp(-cs*f*t_prop)
    return cg_amp

# |G(rn, f)|^2: Practical use
def green_function_ft_mod2(f, cs, t_prop):
    cg_mod2=np.exp(-2.*cs*f*t_prop)
    return cg_mod2

# ---------------------------- #
# ---------- Weight ---------- #
# ---------------------------- #
def weight_calc(r, gamma, r_in, stress):
    #Stressed
    if stress==1:
        w=r**(-gamma)

    #Stress-free
    elif stress==2:
        w=(1.-np.sqrt(r_in/r))*(r**(-gamma))*(r>r_in)
    w/=np.sum(w) # Normalize so that \sum w=1
    return w

def weight_calc_pivot(r, cc, gamma, r_in, stress):
    alpha=gamma-2 # r_n * dr_n \propto r^{+2}
    #Stressed
    if stress==1:
        w=cc*r**(-alpha)

    #Stress-free
    elif stress==2:
        w=cc*(1.-np.sqrt(r_in/r))*(r**(-alpha))*(r>r_in)
    #w/=np.sum(w) # Normalize so that \sum w=1
    return w
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

# Normalized transfer function for reprocessing
def ch_rep_calc(f, t0, dt0):
    omega=2.*np.pi*f
    lh=np.sinc(omega*dt0/(2.*np.pi))*np.exp(-1j*omega*t0)
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

# No reverberation
def lf2_calc(fs,\
             n_r,\
             ws,\
             lm2s,\
             fs_prop,\
             cs,\
             dr_r,\
             rg_c):
    ts_prop=1./fs_prop #[Rg/c]
    tot=0
    for i_r in range(n_r):
        tot+=(ws[i_r]**2)*lm2s[i_r]

        if i_r==0:
            continue

        # Cross term
        tot_c=0
        for i_ro in range(i_r):
            ### Propagation time ###
            t_prop=prop_time_calc(i_start=i_ro, i_end=i_r, ts_vis=ts_prop, dr_r=dr_r) #[Rg/c]
            t_prop*=rg_c #[s]

            # |G(r_k, r_n, f)|
            cgs_amp=green_function_ft_amp(f=fs, cs=cs, t_prop=t_prop)

            ### Cross term ###
            tot_c+=ws[i_ro]*ws[i_r]*np.cos(2.*np.pi*fs*t_prop)*lm2s[i_ro]*cgs_amp

        tot+=2.*tot_c

    return tot

# Reverberation included
def lf2_ref_calc(fs,\
                 n_r,\
                 ws,\
                 wrss,\
                 lm2s,\
                 fs_prop,\
                 cs,\
                 dr_r,\
                 rg_c):
    # Propagation time
    ts_prop=1./fs_prop #[Rg/c]

    # Weight
    wrs_tmp=wrss[0]
    wrss_tmp, wss=np.meshgrid(wrs_tmp, ws)
    #print(wss)
    #print(wss.shape)
    wss=wss+wrss
    #print(wss)
    #print(wss.shape)

    tot=0
    for i_r in range(n_r):
        tot+=(np.abs(wss[i_r])**2)*lm2s[i_r]

        if i_r==0:
            continue

        # Cross term
        tot_c=0
        for i_ro in range(i_r):
            ### Propagation time ###
            t_prop=prop_time_calc(i_start=i_ro, i_end=i_r, ts_vis=ts_prop, dr_r=dr_r) #[Rg/c]
            t_prop*=rg_c #[s]

            # |G(r_k, r_n, f)|
            cgs_amp=green_function_ft_amp(f=fs, cs=cs, t_prop=t_prop)

            ### Cross term ###
            tot_c+=np.real(np.conjugate(wss[i_ro])*wss[i_r]*np.exp(-1j*2.*np.pi*fs*t_prop))*lm2s[i_ro]*cgs_amp

        tot+=2.*tot_c

    return tot

# No reverberation
def lflf_calc(fs,\
              n_r,\
              ws,\
              ws_r,\
              lm2s,\
              fs_prop,\
              cs,\
              dr_r,\
              rg_c):
    ts_prop=1./fs_prop #[Rg/c]
    tot=0
    for i_r in range(n_r):
        tot+=ws_r[i_r]*ws[i_r]*lm2s[i_r]

        if i_r==0:
            continue

        # Cross term
        tot_c=0
        for i_ro in range(i_r):
            ### Propagation time ###
            t_prop=prop_time_calc(i_start=i_ro, i_end=i_r, ts_vis=ts_prop, dr_r=dr_r) #[Rg/c]
            t_prop*=rg_c #[s]

            # |G(r_k, r_n, f)|
            cgs_amp=green_function_ft_amp(f=fs, cs=cs, t_prop=t_prop)

            ### Cross term ###
            tot_c+=(ws_r[i_ro]*ws[i_r] *np.exp(-1j*2.*np.pi*fs*t_prop) + \
                    ws_r[i_r] *ws[i_ro]*np.exp(1j*2.*np.pi*fs*t_prop ))*lm2s[i_ro]*cgs_amp

        tot=tot+tot_c

    return tot

# Reverberation included
def lflf_ref_calc(fs,\
                  n_r,\
                  ws,\
                  wrss,\
                  ws_r,\
                  wrss_r,\
                  lm2s,\
                  fs_prop,\
                  cs,\
                  dr_r,\
                  rg_c):
    # Propagation time
    ts_prop=1./fs_prop #[Rg/c]

    # Weight
    # --- Subject band --- #
    wrs_tmp=wrss[0]
    wrss_tmp, wss=np.meshgrid(wrs_tmp, ws)
    wss=wss+wrss
    # --- Reference band --- #
    wrs_r_tmp=wrss_r[0]
    wrss_r_tmp, wss_r=np.meshgrid(wrs_r_tmp, ws_r)
    wss_r=wss_r+wrss_r

    tot=0
    for i_r in range(n_r):
        tot+=np.conjugate(wss_r[i_r])*wss[i_r]*lm2s[i_r]

        if i_r==0:
            continue

        # Cross term
        tot_c=0
        for i_ro in range(i_r):
            ### Propagation time ###
            t_prop=prop_time_calc(i_start=i_ro, i_end=i_r, ts_vis=ts_prop, dr_r=dr_r) #[Rg/c]
            t_prop*=rg_c #[s]

            # |G(r_k, r_n, f)|
            cgs_amp=green_function_ft_amp(f=fs, cs=cs, t_prop=t_prop)

            ### Cross term ###
            tot_c+=(np.conjugate(wss_r[i_ro])*wss[i_r] *np.exp(-1j*2.*np.pi*fs*t_prop) + \
                    np.conjugate(wss_r[i_r]) *wss[i_ro]*np.exp(1j*2.*np.pi*fs*t_prop ))*lm2s[i_ro]*cgs_amp

        tot=tot+tot_c

    return tot

def print_ring_info(name, xs, digit):
    print('--------------------------------------------------')
    print('{0}: '.format(name), end='')
    n_x=len(xs)
    if n_x==0:
        print('')
    else:
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

