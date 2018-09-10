def atoms_inspheres(rs, stloc):

    atom_list =['H','C','S','Ag','Au']

    rs_pos = stloc.positions[rs]
    dists = dist_between(p1 = stloc.positions,p2 = rs_pos)
    stsym = np.array(stloc.symbols)

    bondsl = []
    atcount = []
    num_layers = 4
    for i in range(0,num_layers):
        atomloc = [np.where(np.logical_and(dists>= i*3.5 , dists< (i+1)*3.5))[0]] 
        g = stsym[atomloc]
        apos = stloc.positions[atomloc]

        Hn =  np.count_nonzero(g == 'H')
        Cn = np.count_nonzero(g == 'C')
        Sn = np.count_nonzero(g == 'S')
        Agn = np.count_nonzero(g == 'Ag')
        Aun = np.count_nonzero(g == 'Au')
    #     print(Hn+Sn+Cn+Agn+Aun)
    #     print(Hn,Sn,Cn,Agn,Aun)
    #   Number of atoms in different spheres
        atcount.append( [Hn,Sn,Cn,Agn,Aun, Hn+Sn+Cn+Agn+Aun] )

        pos_a =[]
        for j in range(len(atom_list)):
            al = atom_list[j]
            curr_a = apos[g == al]
            pos_a.append(curr_a)
    #     atom_list =['H','C','S','Ag','Au']
    #                 0    1    2   3    4

        # H-C bonds = 1.14 +/- .2
        HCd = 1.14
        HCbonds = 0
        if len(pos_a[1]):
            for il in range(len(pos_a[1])):
                d = dist_between(pos_a[0], pos_a[1][il]) 
                c = np.count_nonzero(np.logical_and(d < HCd +.2, d > HCd -.2 ))
                HCbonds += c

    #     C-S bonds = 1.9 +/- .2
        CSd = 1.9
        CSbonds = 0
        if len(pos_a[2]):
            for il in range( len(pos_a[2]) ):
                d = dist_between(pos_a[1],pos_a[2][il]) 
                c = np.count_nonzero(np.logical_and(d < CSd +.2, d > CSd -.2 ))
                CSbonds += c

        # S-Ag bonds = 2.4 +/- .2
        SAgd = 2.4
        SAgbonds = 0
        if len(pos_a[3]):
            for il in range(len(pos_a[3])):
                d = dist_between(pos_a[2], pos_a[3][il]) 
                c = np.count_nonzero(np.logical_and(d < SAgd +.2, d > SAgd -.2 ))
                SAgbonds += c

        # S-Au bonds = 2.4 +/- .2
        SAud = 2.4
        SAubonds = 0
        if len(pos_a[4]):
            for il in range(len(pos_a[4])):
                d = dist_between(pos_a[2], pos_a[4][il]) 
                c = np.count_nonzero(np.logical_and(d < SAud +.2, d > SAud -.2 ))
                SAubonds += c
    #     print(SAubonds)
    #     atom_list =['H','C','S','Ag','Au']
    #                 0    1    2   3    4
        # Ag-Ag bonds = 3.5 +/- .5
        AgAgd = 3.5
        AgAgbonds = 0
        if len(pos_a[3]):
            for il in range(len(pos_a[3])):
                d = dist_between(pos_a[3], pos_a[3][il]) 
                c = np.count_nonzero(np.logical_and(d < AgAgd +.5, d > AgAgd -.5 ))
                AgAgbonds += c

        # Ag-Au bonds = 2.9 +/- .2
        AgAud = 2.9
        AgAubonds = 0
        if len(pos_a[3]):
            for il in range(len(pos_a[3])):
                d = dist_between(pos_a[4], pos_a[3][il]) 
                c = np.count_nonzero(np.logical_and(d < AgAud +.2, d > AgAud -.2 ))
                AgAubonds += c

        # Au-Au bonds  = 2.8 +/- .2
        AuAud = 2.9
        AuAubonds = 0
        if len(pos_a[4]):
            for il in range(len(pos_a[4])):
                d = dist_between(pos_a[4], pos_a[4][il]) 
                c = np.count_nonzero(np.logical_and(d < AuAud +.2, d > AuAud -.2 ))
                AuAubonds += c
        bondsl.append( [HCbonds, CSbonds,SAgbonds, SAubonds, AgAgbonds, AgAubonds, AuAubonds ])
    # #         print(bondsl)
    #         allbondsl.append(bondsl)

    bl = [bondsl[ix] for ix in range(num_layers)]
    bl = np.array(bl)
    bl = bl.reshape(1,bl.shape[0]*bl.shape[1])

    atc = [atcount[ix] for ix in range(num_layers)]
    atc = np.array(atc)
    atc = atc.reshape(1,atc.shape[0]*atc.shape[1])

    return([atc, bl])

