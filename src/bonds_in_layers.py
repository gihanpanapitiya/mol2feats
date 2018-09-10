def B1_B2count10(p_lc , s_lc, startl, endl):

    HCd = 1.14
    CSd = 1.9
    SAgd = 2.4
    SAud = 2.4
    AgAgd = 3.5
    AgAud = 2.9
    AuAud = 2.9

    lp = p_lc[startl:endl]

    ls = s_lc[startl:endl]

    #HC
    Hpos = np.where(ls == 1.)[0]
    Cpos = np.where(ls == 6.)[0]
    Spos = np.where(ls == 16.)[0]
    Agpos = np.where(ls == 47.)[0]
    Aupos = np.where(ls == 79.)[0]

    HCbonds = 0
    if Cpos.size:
        for il in range(len(Cpos)):
            d = dist_between(lp[Hpos], lp[Cpos[il]]) 
            c = np.count_nonzero(np.logical_and(d < HCd +.2, d > HCd -.2 ))
            HCbonds += c

    CSbonds = 0
    if Spos.size:
        for il in range(len(Spos)):
            d = dist_between(lp[Cpos], lp[Spos[il]]) 
            c = np.count_nonzero(np.logical_and(d < CSd +.2, d > CSd -.2 ))
            CSbonds += c

    SAgbonds = 0
    if Agpos.size:
        for il in range(len(Agpos)):
            d = dist_between(lp[Spos], lp[Agpos[il]]) 
            c = np.count_nonzero(np.logical_and(d < SAgd +.2, d > SAgd -.2 ))
            SAgbonds += c

    SAubonds = 0
    if Aupos.size:
        for il in range(len(Aupos)):
            d = dist_between(lp[Spos], lp[Aupos[il]]) 
            c = np.count_nonzero(np.logical_and(d < SAud +.2, d > SAud -.2 ))
            SAubonds += c

    AgAgbonds = 0
    if Agpos.size:
        for il in range(len(Agpos)):
            d = dist_between(lp[Agpos], lp[Agpos[il]]) 
            c = np.count_nonzero(np.logical_and(d < AgAgd +.5, d > AgAgd -.5 ))
            AgAgbonds += c

    AgAubonds = 0
    if Agpos.size:
        for il in range(len(Agpos)):
            d = dist_between(lp[Aupos], lp[Agpos[il]]) 
            c = np.count_nonzero(np.logical_and(d < AgAud +.2, d > AgAud -.2 ))
            AgAubonds += c

    AuAubonds = 0
    if Aupos.size:
        for il in range(len(Aupos)):
            d = dist_between(lp[Aupos], lp[Aupos[il]]) 
            c = np.count_nonzero(np.logical_and(d < AuAud +.2, d > AuAud -.2 ))
            AuAubonds += c

    bsum = HCbonds + CSbonds + SAgbonds + SAubonds + AgAgbonds + AgAubonds + AuAubonds
    if bsum != 0:
        norm_bsum = [HCbonds/bsum , CSbonds/bsum , SAgbonds/bsum , SAubonds/bsum ,\
         AgAgbonds/bsum , AgAubonds/bsum , AuAubonds/bsum]
    else:
        norm_bsum = [0 , 0 , 0 , 0 , 0 , 0 , 0]

    return(norm_bsum)
