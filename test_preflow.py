from preflow import *

def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    engs=np.arange(1.e-2, 1.e1, 1.e-2)

    kTbbd    =0.2
    normd    =5.e4
    Gammas   =1.8
    kTbbs    =kTbbd
    kTes     =1.e2
    norms    =2.
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
    Frefs    =0.01
    Rinh     =45.
    Routh    =400.
    Indexh   =3.
    logxih   =3.
    logNh    =15.
    Frefh    =0.01
    Mass     =8.
    rin      =6.
    drhc     =10.
    drsc     =16.
    drd      =0.
    Nring    =40
    Fvdisk   =0.8
    drdisk   =1.e3
    Fvflow   =0.8
    drflow   =1.e3
    Bdisk    =0.05
    mdisk    =0.3
    Badisk   =0.05
    madisk   =0.3
    Bflow    =4.
    mflow    =1.
    Baflow   =4.
    maflow   =1.
    Sm       =1.e-1
    gammad   =3.
    gammaf   =3.
    stress   =2.
    rmin     =6.
    Emin     =0.5
    Emax     =1.
    #Emin     =1.
    #Emax     =2.6
    Eminr    =2.6
    Emaxr    =4.8
    eta0d    =1.
    eta1d    =0.2
    eta0sc   =1.
    eta1sc   =0.3
    eta0hc   =1.
    eta1hc   =0.4
    eta0sr   =1.
    eta1sr   =-0.2
    eta0hr   =1.
    eta1hr   =-0.5
    t0sc     =1.e-3
    dt0sc    =1.e-2
    t0hc     =1.e-3
    dt0hc    =1.e-2
    quant    =5
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
        Frefs    ,\
        Rinh     ,\
        Routh    ,\
        Indexh   ,\
        logxih   ,\
        logNh    ,\
        Frefh    ,\
        Mass     ,\
        rin      ,\
        drhc     ,\
        drsc     ,\
        drd      ,\
        Nring    ,\
        Fvdisk   ,\
        drdisk   ,\
        Fvflow   ,\
        drflow   ,\
        Bdisk    ,\
        mdisk    ,\
        Badisk   ,\
        madisk   ,\
        Bflow    ,\
        mflow    ,\
        Baflow   ,\
        maflow   ,\
        Sm       ,\
        gammad   ,\
        gammaf   ,\
        stress   ,\
        rmin     ,\
        Emin     ,\
        Emax     ,\
        Eminr    ,\
        Emaxr    ,\
        eta0d    ,\
        eta1d    ,\
        eta0sc   ,\
        eta1sc   ,\
        eta0hc   ,\
        eta1hc   ,\
        eta0sr   ,\
        eta1sr   ,\
        eta0hr   ,\
        eta1hr   ,\
        t0sc     ,\
        dt0sc    ,\
        t0hc     ,\
        dt0hc    ,\
        quant    ,\
        invert   ,\
        par_print\
        ]

    fluxes=np.ones(len(engs)-1)

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    preflow(engs=engs, params=params, fluxes=fluxes)

if __name__=='__main__':
    main()
