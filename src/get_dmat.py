def get_dmat(stf):
    dmatrix = np.zeros((natoms,natoms))

    for j in range(natoms):

        d = dist_between(p1=stf.positions, p2=stf.positions[j])
        s = np.array(stf.symbols)

        for i in range(0,natoms):
            dideal = dict_ari[s[j]] + dict_ari[ s[i] ] +.25
            if (d[i] < dideal and d[i]>0):
                dmatrix[j, i ] = 1
                
    return dmatrix
