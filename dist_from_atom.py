def distfromrs(reactive_site, structure):
# ############################
    #     For Ag and Au
# ############################

#     mol1 = np.loadtxt(str(fpath)+'answer.bas', skiprows=1)

#     mol1_pos = mol1[:,1:4]
#     mol1_atms = mol1[:,0]

#     st_local = pychemia.Structure(positions=mol1_pos, symbols=mol1_atms, periodicity=False)
    st_local = structure
    sa = pychemia.analysis.StructureAnalysis(st_local)

    alldist = sa.all_distances()
    st_syms = np.array(st.symbols)
    auagsites = np.where(np.logical_or(st_syms == 'Au',st_syms == 'Ag'))
    agsites = np.where(st_syms=='Ag')
    

    cv = int(reactive_site - 1)  # Python starts with zero
    dist_list = []
    for i in range(len(auagsites[0])):
        cs = auagsites[0][i]
        if cs < cv:
    #         print (alldist[(cs,cv)], cs, cv)
            dist_list.append([alldist[(cs,cv)], cs, cv])
        else:
    #         print (alldist[(cv,cs)], cs, cv)
            dist_list.append([alldist[(cv,cs)], cs, cv])

    dist = []
    dist_pos = []
    for i in range(len(dist_list)):
        dist.append(dist_list[i][0])
        dist_pos.append(dist_list[i][1])
        
    args = np.argsort(dist) 
    agsites = np.array(agsites)
    asc_sites = np.array(dist_pos)[args]
#     print(asc_sites)
    asc_dists = np.array(dist)[args]
    
#     xx = np.nonzero(np.in1d( asc_sites,agsites ))[0]
#     print(asc_dists[1:len(asc_dists)])
    dlist = [cv+1,asc_dists[1:len(asc_dists)]]
    return dlist
