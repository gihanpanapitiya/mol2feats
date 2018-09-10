def threeBcount(p_lc , s_lc, startl, endl):    
    
    lp = p_lc[startl:endl]
    ls = s_lc[startl:endl]
    
    
    #HC
    Hpos = np.where(ls == 1.)[0]
    Cpos = np.where(ls == 6.)[0]
    Spos = np.where(ls == 16.)[0]
    Agpos = np.where(ls == 47.)[0]
    Aupos = np.where(ls == 79.)[0]
    
    SCH_bonds, CSAu_bonds, AgAuAu_bonds, AuAgAu_bonds, HCH_bonds= 0,0,0,0,0

    SCH_bonds = three_ABC(l=1, c=6, r=16,  lpos=Hpos, cpos=Cpos, rpos=Spos, lp=lp)
    CSAu_bonds = three_ABC(l=6, c=16, r=79,  lpos=Cpos, cpos=Spos, rpos=Aupos, lp=lp)
    AgAuAu_bonds = three_ABC(l=47, c=79, r=79,  lpos=Agpos, cpos=Aupos, rpos=Aupos, lp=lp)
    AuAgAu_bonds = three_ABA(l=79, c=47, r=79,  lpos=Aupos, cpos=Agpos, rpos=Aupos, lp=lp)
    HCH_bonds = three_ABA(l=1, c=6, r=1,  lpos=Hpos, cpos=Cpos, rpos=Hpos, lp=lp)
    
    bsum = SCH_bonds + CSAu_bonds+ AgAuAu_bonds+AuAgAu_bonds+HCH_bonds
#     print(SCH_bonds, CSAu_bonds, AgAuAu_bonds, AuAgAu_bonds, HCH_bonds)
#     print(bsum)
    if bsum != 0:
        x = [ SCH_bonds/bsum , CSAu_bonds/bsum,  AgAuAu_bonds/bsum, AuAgAu_bonds/bsum, HCH_bonds/bsum ]
    else:
        x = [0,0,0,0,0]
    return(x)
