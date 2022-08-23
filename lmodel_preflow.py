##########################################################################
##### Type 'AllModels.addPyMod(preflow, lmodel_preflow_info, 'add')' #####
##### for XSPEC to read the model, preflow.                          #####
##########################################################################
lmodel_preflow_info=\
('Mass    Msun   8.    1.    1.     100.   100.   -1     ',\
 'rin     Rg     6.    1.    1.     1.e2   1.e2   -0.1   ',\
 'drh     Rg     10.   0.    0.     1.e2   1.e2   -0.1   ',\
 'drs     Rg     8.    0.    0.     1.e2   1.e2   -0.1   ',\
 'drd     Rg     13.   0.    0.     1.e3   1.e3   -0.1   ',\
 'Nring   \'\'   40    1     1      100    100    -1     ',\
 'Fvd     \'\'   0.8   0.    0.     1.e1   1.e1   0.01   ',\
 'drvd    \'\'   1.e3  0.1   0.1    1.e3   1.e3   -1.    ',\
 'Fvf     \'\'   0.8   0.    0.     1.e1   1.e1   0.01   ',\
 'drvf    \'\'   1.e3  0.1   0.1    1.e3   1.e3   -1.    ',\
 'Bd      \'\'   0.03  1.e-2 1.e-2  1.e3   1.e3   -1.e-2 ',\
 'md      \'\'   0.5   0.1   0.1    10.    10.    -1.e-2 ',\
 'Bpd     \'\'   0.03  1.e-2 1.e-2  1.e3   1.e3   -1.e-2 ',\
 'mpd     \'\'   0.5   0.1   0.1    10.    10.    -1.e-2 ',\
 'Bf      \'\'   6.    1.e-2 1.e-2  1.e3   1.e3   -1.e-2 ',\
 'mf      \'\'   1.2   0.1   0.1    10.    10.    -1.e-2 ',\
 'Bpf     \'\'   6.    1.e-2 1.e-2  1.e3   1.e3   -1.e-2 ',\
 'mpf     \'\'   1.2   0.1   0.1    10.    10.    -1.e-2 ',\
 'D       \'\'   1.e-1 0.    0.     1.e1   1.e1   -1     ',\
 'indexd  \'\'   3.    0.    0.     10.    10.    -1     ',\
 'indexf  \'\'   3.    0.    0.     10.    10.    -1     ',\
 'stress  \'\'   2     1     1      2      2      -1     ',\
 'Emin    keV    0.5   1.e-3 1.e-3  1.e3   1.e3   -1     ',\
 'Emax    keV    0.6   1.e-3 1.e-3  1.e3   1.e3   -1     ',\
 'Sd      \'\'   0.    -1.   -1.    1.     1.     0.01   ',\
 'Ss      \'\'   0.2   -1.   -1.    1.     1.     0.01   ',\
 'Sh      \'\'   0.2   -1.   -1.    1.     1.     0.01   ',\
 'Ssr     \'\'   0.    -1.   -1.    1.     1.     -1.e-2 ',\
 'Shr     \'\'   0.    -1.   -1.    1.     1.     -1.e-2 ',\
 'Eminr   keV    2.    1.e-3 1.e-3  1.e3   1.e3   -1     ',\
 'Emaxr   keV    10.   1.e-3 1.e-3  1.e3   1.e3   -1     ',\
 'Srd     \'\'   0.    -1.   -1.    1.     1.     0.01   ',\
 'Srs     \'\'   0.2   -1.   -1.    1.     1.     0.01   ',\
 'Srh     \'\'   0.2   -1.   -1.    1.     1.     0.01   ',\
 'Srsr    \'\'   0.    -1.   -1.    1.     1.     -1.e-2 ',\
 'Srhr    \'\'   0.    -1.   -1.    1.     1.     -1.e-2 ',\
 't0s     s      5.e-3 1.e-4 1.e-4  1.     1.     -1     ',\
 'dt0s    s      5.e-3 1.e-4 1.e-4  1.     1.     -1     ',\
 't0h     s      5.e-3 1.e-4 1.e-4  1.     1.     -1     ',\
 'dt0h    s      5.e-3 1.e-4 1.e-4  1.     1.     -1     ',\
 'quant   \'\'   1     1     1      6      6      -1     ',\
 'invert  \'\'   2     1     1      2      2      -1     ',\
 'print   \'\'   1     1     1      2      2      -1     ')

################################################################################
##### Type 'AllModels.addPyMod(preflowscp, lmodel_preflowscp_info, 'add')' #####
##### for XSPEC to read the model, preflow.                                #####
################################################################################
lmodel_preflowscp_info=\
('kTbbd    keV    0.2   1.e-2  1.e-2  1.e2   1.e2   1.e-2  ',\
 'normd    \'\'   3.e5  0.     0.     1.e8   1.e8   1.e4   ',\
 'Gammas   \'\'   2.    1.001  1.001  10.    10.    1.e-2  ',\
 'kTbbs    keV    0.2   1.e-2  1.e-2  1.e1   1.e1   1.e-2  ',\
 'kTes     keV    100.  1.e0   1.e0   1.e3   1.e3   1.e-1  ',\
 'norms    \'\'   2.    0.     0.     1.e2   1.e2   1.e-2  ',\
 'Gammah   \'\'   1.5   1.001  1.001  10.    10.    1.e-2  ',\
 'kTbbh    keV    0.2   1.e-2  1.e-2  1.e1   1.e1   1.e-2  ',\
 'kTeh     keV    100.  1.e0   1.e0   1.e3   1.e3   1.e-1  ',\
 'normh    \'\'   2.    0.     0.     1.e2   1.e2   1.e-2  ',\
 'Incl     deg    60.   0.     0.     90.    90.    -1.e0  ',\
 'a        \'\'   0.    -0.998 -0.998 0.998  0.998  -1.e-2 ',\
 'Afe      \'\'   1.    0.5    0.5    10.    10.    -1.e-1 ',\
 'Rins     Rg     60.   1.     1.     1.e3   1.e3   -1.e0  ',\
 'Routs    Rg     4.e2  1.     1.     1.e3   1.e3   -1.e0  ',\
 'Indexs   \'\'   3.    0.     0.     10.    10.    -1.e-1 ',\
 'logxis   \'\'   2.5   0.     0.     4.7    4.7    1.e-2  ',\
 'logNs    cm^-3  15.   15.    15.    20.    20.    -1.e0  ',\
 'normsr   \'\'   0.    0.     0.     1.e1   1.e1   1.e-2  ',\
 'Rinh     Rg     60.   1.     1.     1.e3   1.e3   -1.e0  ',\
 'Routh    Rg     4.e2  1.     1.     1.e3   1.e3   -1.e0  ',\
 'Indexh   \'\'   3.    0.     0.     10.    10.    -1.e-1 ',\
 'logxih   \'\'   2.5   0.     0.     4.7    4.7    1.e-2  ',\
 'logNh    cm^-3  15.   15.    15.    20.    20.    -1.e0  ',\
 'normhr   \'\'   0.    0.     0.     1.e1   1.e1   1.e-2  ',\
 'Mass     Msun   8.    1.     1.     100.   100.   -1.e-1 ',\
 'rin      Rg     6.    1.     1.     1.e2   1.e2   -0.1   ',\
 'drh      Rg     10.   0.     0.     1.e2   1.e2   -0.1   ',\
 'drs      Rg     16.   0.     0.     1.e2   1.e2   -0.1   ',\
 'drd      Rg     13.   0.     0.     1.e3   1.e3   -0.1   ',\
 'Nring    \'\'   40    1      1      100    100    -1     ',\
 'Fvd      \'\'   0.8   0.     0.     1.e1   1.e1   0.01   ',\
 'drvd     \'\'   1.e3  0.1    0.1    1.e3   1.e3   -1.    ',\
 'Fvf      \'\'   0.8   0.     0.     1.e1   1.e1   0.01   ',\
 'drvf     \'\'   1.e3  0.1    0.1    1.e3   1.e3   -1.    ',\
 'Bd       \'\'   0.03  1.e-2  1.e-2  1.e3   1.e3   -1.e-2 ',\
 'md       \'\'   0.5   0.1    0.1    10.    10.    -1.e-2 ',\
 'Bpd      \'\'   0.03  1.e-2  1.e-2  1.e3   1.e3   -1.e-2 ',\
 'mpd      \'\'   0.5   0.1    0.1    10.    10.    -1.e-2 ',\
 'Bf       \'\'   4.    1.e-2  1.e-2  1.e3   1.e3   -1.e-2 ',\
 'mf       \'\'   1.    0.1    0.1    10.    10.    -1.e-2 ',\
 'Bpf      \'\'   4.    1.e-2  1.e-2  1.e3   1.e3   -1.e-2 ',\
 'mpf      \'\'   1.    0.1    0.1    10.    10.    -1.e-2 ',\
 'D        \'\'   0.    0.     0.     1.e1   1.e1   -1     ',\
 'indexd   \'\'   3.    0.     0.     10.    10.    -1     ',\
 'indexf   \'\'   3.    0.     0.     10.    10.    -1     ',\
 'stress   \'\'   2     1      1      2      2      -1     ',\
 'Emin     keV    35.   1.e-3  1.e-3  1.e3   1.e3   -1     ',\
 'Emax     keV    48.   1.e-3  1.e-3  1.e3   1.e3   -1     ',\
 'Eminr    keV    2.4   1.e-3  1.e-3  1.e3   1.e3   -1     ',\
 'Emaxr    keV    4.8   1.e-3  1.e-3  1.e3   1.e3   -1     ',\
 'eta0d    \'\'   1.    -1.    -1.    1.     1.     1.e-2  ',\
 'eta1d    \'\'   0.    -1.    -1.    1.     1.     1.e-2  ',\
 'eta0s    \'\'   1.    -1.    -1.    1.     1.     1.e-2  ',\
 'eta1s    \'\'   0.    -1.    -1.    1.     1.     1.e-2  ',\
 'eta0h    \'\'   1.    -1.    -1.    1.     1.     1.e-2  ',\
 'eta1h    \'\'   0.    -1.    -1.    1.     1.     1.e-2  ',\
 'eta0sr   \'\'   1.    -1.    -1.    1.     1.     1.e-2  ',\
 'eta1sr   \'\'   0.    -1.    -1.    1.     1.     1.e-2  ',\
 'eta0hr   \'\'   1.    -1.    -1.    1.     1.     1.e-2  ',\
 'eta1hr   \'\'   0.    -1.    -1.    1.     1.     1.e-2  ',\
 't0s      s      1.e-3 1.e-4  1.e-4  1.     1.     -1     ',\
 'dt0s     s      1.e-2 1.e-4  1.e-4  1.     1.     -1     ',\
 't0h      s      1.e-3 1.e-4  1.e-4  1.     1.     -1     ',\
 'dt0h     s      1.e-2 1.e-4  1.e-4  1.     1.     -1     ',\
 'quant    \'\'   1     0      0      6      6      -1     ',\
 'invert   \'\'   2     1      1      2      2      -1     ',\
 'print    \'\'   1     1      1      2      2      -1     ')

