from preflow import *

def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    engs=np.arange(1.e-2, 1.e1, 1.e-2)

    Mass  =8.
    rin   =6.
    drh   =10.
    drs   =16.
    drd   =13.
    Nring =40.
    Fvd   =0.8
    drvd  =1.e3
    Fvf   =0.8
    drvf  =1.e3
    Bd    =0.03
    md    =0.5
    Bpd   =0.03
    mpd   =0.5
    Bf    =4.
    mf    =1.
    Bpf   =4.
    mpf   =1.
    D     =0.
    indexd=3.
    indexf=3.
    stress=2.
    Emin  =0.5
    Emax  =1.
    Sd    =0.1
    Ss    =0.2
    Sh    =0.2
    Ssr   =0.1
    Shr   =0.1
    Sps   =0.2
    Sph   =0.2
    Spsr  =0.1
    Sphr  =0.1
    Eminr =2.6
    Emaxr =4.8
    Srd   =0.
    Srs   =0.4
    Srh   =0.4
    Srsr  =0.
    Srhr  =0.
    Srps  =0.4
    Srph  =0.4
    Srpsr =0.
    Srphr =0.
    t0s   =1.e-3
    dt0s  =1.e-2
    t0h   =1.e-3
    dt0h  =1.e-2
    quant =1
    invert=2
    par_print=1

    params=[\
        Mass  ,\
        rin   ,\
        drh   ,\
        drs   ,\
        drd   ,\
        Nring ,\
        Fvd   ,\
        drvd  ,\
        Fvf   ,\
        drvf  ,\
        Bd    ,\
        md    ,\
        Bpd   ,\
        mpd   ,\
        Bf    ,\
        mf    ,\
        Bpf   ,\
        mpf   ,\
        D     ,\
        indexd,\
        indexf,\
        stress,\
        Emin  ,\
        Emax  ,\
        Sd    ,\
        Ss    ,\
        Sh    ,\
        Ssr   ,\
        Shr   ,\
        Srs   ,\
        Srh   ,\
        Srsr  ,\
        Srhr  ,\
        Eminr ,\
        Emaxr ,\
        Srd   ,\
        Srs   ,\
        Srh   ,\
        Srsr  ,\
        Srhr  ,\
        Srps  ,\
        Srph  ,\
        Srpsr ,\
        Srphr ,\
        t0s   ,\
        dt0s  ,\
        t0h   ,\
        dt0h  ,\
        quant ,\
        invert,\
        par_print,\
        ]

    fluxes=np.ones(len(engs)-1)

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    preflow(engs=engs, params=params, fluxes=fluxes)

if __name__=='__main__':
    main()

