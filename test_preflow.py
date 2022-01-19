from preflow import *

def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    engs=np.arange(1.e-2, 1.e1, 1.e-2)

    Mass   =8.
    rin    =6.
    rmh    =16.
    rsm    =24.
    rds    =32.
    rout   =45.
    Nring  =40
    tref   =4.e-3
    dtref  =4.e-3
    Fvdisk =0.8
    drdisk =100.
    Fvflow =0.8
    drflow =100.
    Bdisk  =0.03
    mdisk  =0.5
    Bflow  =6.
    mflow  =1.2
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
    Cd     =0.
    Csc    =1.
    Cmc    =1.
    Chc    =1.
    Eminr  =2.
    Emaxr  =10.
    Cdr    =0.
    Cscr   =-1.
    Cmcr   =1.
    Chcr   =1.
    Eminrr =0.5
    Emaxrr =10.
    Cdrr   =0.
    Cscrr  =1.
    Cmcrr  =-1.
    Chcrr  =1.
    quant  =2
    invert =2
    par_print=1 

    params=[\
        Mass,   rin,    rmh,    rsm,    rds,\
        rout,   Nring,  tref,   dtref,  Fvdisk,\
        drdisk, Fvflow, drflow, Bdisk,  mdisk,\
        Bflow,  mflow,  Ddisk,  Dflow,  Dtran,\
        xlag,   gammad, gammaf, stress, rmin,\
        Emin,   Emax,   Cd,     Csc,    Cmc,\
        Chc,    Eminr,  Emaxr,  Cdr,    Cscr,\
        Cmcr,   Chcr,   Eminrr, Emaxrr, Cdrr,\
        Cscrr,  Cmcrr,  Chcrr,  quant,  invert,\
        par_print]

    fluxes=np.ones(len(engs)-1)

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    preflow(engs=engs, params=params, fluxes=fluxes)

if __name__=='__main__':
    main()
