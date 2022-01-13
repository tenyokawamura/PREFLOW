from preflow import *

def main():
    # ----------------------------- #
    # ---------- Setting ---------- #
    # ----------------------------- #
    engs=np.arange(1.e-2, 1.e1, 1.e-2)
    #engs=np.arange(1.e-2, 1.e1, 1.e-3)
    # 0.  mass
    # 1.  r_in
    # 2.  r_sh
    # 3.  r_ds
    # 4.  r_out
    # 5.  n_ring
    # 6.  tref
    # 7.  dtref
    # 8.  lf_var
    # 9.  lb_disk
    # 10. m_disk
    # 11. lb_flow
    # 12. m_flow
    # 13. stress
    # 14. gamma
    # 15. E_min
    # 16. E_max
    # 17. frac_disk
    # 18. frac_scomp
    # 19. frac_hcomp
    # 20. frac_sref
    # 21. frac_href
    # 22. E_minr
    # 23. E_maxr
    # 24. frac_diskr
    # 25. frac_scompr
    # 26. frac_hcompr
    # 27. frac_srefr
    # 28. frac_hrefr
    # 39. E_minrr
    # 30. E_maxrr
    # 31. frac_scomprr
    # 32. quant
    # 33. display
    params=[8.,\
            6.,\
            24.,\
            32.,\
            45.,\
            40,\
            4.e-3,\
            4.e-3,\
            0.8,\
            100.,\
            0.8,\
            100.,\
            0.03,\
            0.5,\
            6.,\
            1.2,\
            1.,\
            1.,\
            1.,\
            2,\
            0.51,\
            1.50,\
            1.,\
            -1.,\
            1.,\
            6.,\
            2.01,\
            10.,\
            1.,\
            1.,\
            1.,\
            6.,\
            0.51,\
            10.,\
            1.,\
            1.,\
            1.,\
            6.,\
            1,\
            2,\
            1]
    fluxes=np.ones(len(engs)-1)

    # -------------------------- #
    # ---------- Main ---------- #
    # -------------------------- #
    preflow(engs=engs, params=params, fluxes=fluxes)

if __name__=='__main__':
    main()
