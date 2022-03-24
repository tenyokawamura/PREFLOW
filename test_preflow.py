from preflow import *

def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    engs=np.arange(1.e-2, 1.e1, 1.e-2)

    Mass   =8.
    rin    =6.
    drhc   =10.
    drmc   =8.
    drsc   =8.
    drd    =13.
    Nring  =40
    Fvdisk =0.8
    drdisk =100.
    Fvflow =0.8
    drflow =100.
    Bdisk  =0.03
    mdisk  =0.5
    Badisk =0.03
    madisk =0.5
    Bflow  =4.
    mflow  =1.
    Baflow =6.
    maflow =1.2
    Ddisk  =1.
    Dflow  =1.
    Dtran  =1.
    xlag   =1.
    gammad =3.
    gammaf =3.
    stress =2.
    rmin   =6.
    Emin   =0.5
    Emax   =2.
    Cd     =0.1
    Csc    =0.1
    Cmc    =0.25
    Chc    =0.25
    Fd     =0
    Fsc    =0.1
    Fmc    =0.2
    Fhc    =0.3
    Crd    =0.
    Crsc   =0.
    Crmc   =0.
    Crhc   =0.1
    Frd    =0.
    Frsc   =0.
    Frmc   =0.
    Frhc   =0.
    Eminr  =2.
    Emaxr  =10.
    Cdr    =0.
    Cscr   =-1.
    Cmcr   =1.
    Chcr   =1.
    Fdr    =0
    Fscr   =0.2
    Fmcr   =0.4
    Fhcr   =0.6
    Crdr   =0.
    Crscr  =0.
    Crmcr  =0.
    Crhcr  =0.
    Frdr   =0.
    Frscr  =0.
    Frmcr  =0.
    Frhcr  =0.
    trd    =4.e-3
    dt0d   =4.e-3
    trsc   =4.e-3
    dt0sc  =4.e-3
    trmc   =4.e-3
    dt0mc  =4.e-3
    trhc   =4.e-3
    dt0hc  =4.e-3
    quant  =2
    invert =2
    par_print=1 

    params=[\
        Mass,   rin,    drhc,   drmc,   drsc,\
        drd,    Nring,  Fvdisk, drdisk, Fvflow,\
        drflow, Bdisk,  mdisk,  Badisk, madisk,\
        Bflow,  mflow,  Baflow, maflow, Ddisk,\
        Dflow,  Dtran,  xlag,   gammad, gammaf,\
        stress, rmin,   Emin,   Emax,   Cd,\
        Csc,    Cmc,    Chc,    Fd,     Fsc,\
        Fmc,    Fhc,    Crd,    Crsc,   Crmc,\
        Crhc,   Frd,    Frsc,   Frmc,   Frhc,\
        Eminr,  Emaxr,  Cdr,    Cscr,   Cmcr,\
        Chcr,   Fdr,    Fscr,   Fmcr,   Fhcr,\
        Crdr,   Crscr,  Crmcr,  Crhcr,  Frdr,\
        Frscr,  Frmcr,  Frhcr,  trd,    dt0d,\
        trsc,   dt0sc,  trmc,   dt0mc,  trhc,\
        dt0hc,  quant,  invert, par_print]

    fluxes=np.ones(len(engs)-1)

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    preflow(engs=engs, params=params, fluxes=fluxes)

if __name__=='__main__':
    main()
