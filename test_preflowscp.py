from preflowscp import *

def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    engs=np.arange(1.e-2, 1.e1, 1.e-2)

    kTbbd    =0.2
    normd    =0.
    Gammas   =1.8
    kTbbs    =kTbbd
    kTes     =1.e2
    norms    =0.
    Gammah   =1.5
    kTbbh    =kTbbd
    kTeh     =kTes
    normh    =1.5
    Incl     =60.
    a        =0.998
    Afe      =1.
    Rins     =45.
    Routs    =400.
    Indexs   =3.
    logxis   =3.
    logNs    =15.
    normsr   =0.01
    Rinh     =45.
    Routh    =400.
    Indexh   =3.
    logxih   =3.
    logNh    =15.
    normhr   =0.
    Mass     =8.
    rin      =6.
    drh      =26.
    drs      =0.
    drd      =13.
    Nring    =40
    Fvd      =0.8
    Fvs      =0.8
    Fvh      =0.8
    fgdo     =5.e-2
    md       =0.5
    fpdo     =5.e-2
    mpd      =0.5
    fgfo     =1.
    mf       =1.
    fpfo     =1.
    mpf      =1.
    D        =1.e-1
    indexd   =3.
    indexf   =3.
    stress   =2.
    Emin     =35.
    Emax     =48.
    Eminr    =2.6
    Emaxr    =4.8
    eta0d    =1.
    eta10d   =0.2
    eta0s    =1.
    eta10s   =0.3
    eta0h    =1.
    eta10h   =-0.2
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
        kTbbs    ,\
        kTes     ,\
        norms    ,\
        Gammah   ,\
        kTbbh    ,\
        kTeh     ,\
        normh    ,\
        Incl     ,\
        a        ,\
        Afe      ,\
        Rins     ,\
        Routs    ,\
        Indexs   ,\
        logxis   ,\
        logNs    ,\
        normsr   ,\
        Rinh     ,\
        Routh    ,\
        Indexh   ,\
        logxih   ,\
        logNh    ,\
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
        fgdo     ,\
        md       ,\
        fpdo     ,\
        mpd      ,\
        fgfo     ,\
        mf       ,\
        fpfo     ,\
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
        eta10d   ,\
        eta0s    ,\
        eta10s   ,\
        eta0h    ,\
        eta10h   ,\
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
    preflowscp(engs=engs, params=params, fluxes=fluxes)

if __name__=='__main__':
    main()
