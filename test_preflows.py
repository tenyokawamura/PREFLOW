from preflows import *

def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    engs=np.arange(1.e-2, 1.e1, 1.e-2)

    kTbbd    =0.2
    normd    =5.e4
    Gammas   =1.8
    Ecuts    =1.e2
    norms    =2.
    Gammah   =1.5
    Ecuth    =Ecuts
    normh    =1.5
    Incl     =60.
    a        =0.998
    Afe      =1.
    Rins     =45.
    Routs    =400.
    Indexs   =3.
    logxis   =3.
    #normsr   =0.01
    normsr   =0.
    Rinh     =45.
    Routh    =400.
    Indexh   =3.
    logxih   =3.
    #normhr   =0.01
    normhr   =0.
    Mass     =8.
    rin      =6.
    drh      =10.
    drs      =16.
    drd      =13.
    Nring    =40
    Fvd      =0.8
    Fvs      =0.7
    Fvh      =0.6
    Bd       =0.05
    md       =0.3
    Bpd      =0.05
    mpd      =0.3
    Bf       =4.
    mf       =1.
    Bpf      =4.
    mpf      =1.
    D        =1.e-1
    indexd   =3.
    indexf   =3.
    stress   =2.
    Emin     =0.5
    Emax     =1.
    Eminr    =2.6
    Emaxr    =4.8
    eta0d    =1.
    eta1d    =0.2
    eta0s    =1.
    eta1s    =0.3
    eta0h    =1.
    eta1h    =0.4
    etap0s   =1.
    etap1s   =-0.2
    etap0h   =1.
    etap1h   =-0.5
    t0s      =1.e-3
    dt0s     =1.e-2
    t0h      =1.e-3
    dt0h     =1.e-2
    quant    =2
    invert   =2
    par_print=1

    params=[\
        kTbbd    ,\
        normd    ,\
        Gammas   ,\
        Ecuts    ,\
        norms    ,\
        Gammah   ,\
        Ecuth    ,\
        normh    ,\
        Incl     ,\
        a        ,\
        Afe      ,\
        Rins     ,\
        Routs    ,\
        Indexs   ,\
        logxis   ,\
        normsr   ,\
        Rinh     ,\
        Routh    ,\
        Indexh   ,\
        logxih   ,\
        normhr   ,\
        Mass     ,\
        rin      ,\
        drh      ,\
        drs      ,\
        drd      ,\
        Nring    ,\
        Fvd      ,\
        Fvs      ,\
        Fvh      ,\
        Bd       ,\
        md       ,\
        Bpd      ,\
        mpd      ,\
        Bf       ,\
        mf       ,\
        Bpf      ,\
        mpf      ,\
        D        ,\
        indexd   ,\
        indexf   ,\
        stress   ,\
        Emin     ,\
        Emax     ,\
        Eminr    ,\
        Emaxr    ,\
        eta0d    ,\
        eta1d    ,\
        eta0s    ,\
        eta1s    ,\
        eta0h    ,\
        eta1h    ,\
        etap0s   ,\
        etap1s   ,\
        etap0h   ,\
        etap1h   ,\
        t0s      ,\
        dt0s     ,\
        t0h      ,\
        dt0h     ,\
        quant    ,\
        invert   ,\
        par_print\
        ]

    fluxes=np.ones(len(engs)-1)

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    preflows(engs=engs, params=params, fluxes=fluxes)

if __name__=='__main__':
    main()
