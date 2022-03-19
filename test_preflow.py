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
    Fdd    =0
    Fdsc   =0.1
    Fdmc   =0.2
    Fdhc   =0.3
    Arep   =0.
    Eminr  =2.
    Emaxr  =10.
    Cdr    =0.
    Cscr   =-1.
    Cmcr   =1.
    Chcr   =1.
    Fddr   =0
    Fdscr  =0.2
    Fdmcr  =0.4
    Fdhcr  =0.6
    Arepr  =1.
    Eminrr =0.5
    Emaxrr =10.
    Cdrr   =0.
    Cscrr  =1.
    Cmcrr  =0.
    Chcrr  =0.
    tref   =4.e-3
    dtref  =4.e-3
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
        Csc,    Cmc,    Chc,    Fdd,    Fdsc,\
        Fdmc,   Fdhc,   Arep,   Eminr,  Emaxr,\
        Cdr,    Cscr,   Cmcr,   Chcr,   Fddr,\
        Fdscr,  Fdmcr,  Fdhcr,  Arepr,  Eminrr,\
        Emaxrr, Cdrr,   Cscrr,  Cmcrr,  Chcrr,\
        tref,   dtref,  quant,  invert, par_print]

    fluxes=np.ones(len(engs)-1)

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    preflow(engs=engs, params=params, fluxes=fluxes)

if __name__=='__main__':
    main()
