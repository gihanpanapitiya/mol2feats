# 1. Number of Ag atoms up to 3 links of distance = nag_3 
def get_nag3_ntot2(GG, adsite, agposs):
    
    asite = adsite
    pcont = []
    for i in range(natoms):    
        ag_s = i

        sp = nx.shortest_path(G=GG,source=asite,target=ag_s)
        pcont.append(sp)

    pcont = np.array(pcont)

    x = [len(pcont[i]) < 5 for i in range(natoms)] # upto 3 links
    x1 = pcont[np.nonzero(x)[0]]

    x2  = [ list(set(agposs).intersection(set(x1[j])) ) for j in range(len(x1))]

    flat_list = [item for sublist in x2 for item in sublist]

    nag_3 = len(np.unique(flat_list))
    
# 2. Total number of atoms up to 2 links of distance = ntot_2
#     For ntot_2
    x = [len(pcont[i]) < 4 for i in range(natoms)] # upto 2 links
    x1 = pcont[np.nonzero(x)[0]]

    # x2  = [ list(set(agsites[0]).intersection(set(x1[j])) ) for j in range(len(x1))]

    flat_list = [item for sublist in x1 for item in sublist]

    ntot_2 = len(np.unique(flat_list))

    
    
    
    return [nag_3, ntot_2]
