import networkx as nx
import numpy as np
import itertools


# 1. Number of Ag atoms up to 3 links of distance = nag_3
def get_nlinks_atom(GG, site, atom_locs, upto):
    natoms = len(GG)

    asite = site
    pcont = []
    for i in range(natoms):
        ag_s = i

        sp = nx.shortest_path(G = GG, source = asite, target = ag_s)
        pcont.append(sp)

    pcont = np.array(pcont)

    x = [len(pcont[i]) < upto for i in range(natoms)] # upto 3 links
    x1 = pcont[np.nonzero(x)[0]]

    x2  = [ list(set(atom_locs).intersection(set(x1[j])) ) for j in range(len(x1))]

    flat_list = [item for sublist in x2 for item in sublist]
    nag_3 = len(np.unique(flat_list))


    return nag_3



def get_nlinks_gen(GG, site, upto):
    natoms = len(GG)

    asite = site
    pcont = []
    for i in range(natoms):
        ag_s = i

        sp = nx.shortest_path(G = GG, source = asite, target = ag_s)
        pcont.append(sp)

    pcont = np.array(pcont)

    x = [len(pcont[i]) < upto for i in range(natoms)] # upto 2 links
    x1 = pcont[np.nonzero(x)[0]]

    flat_list = [item for sublist in x1 for item in sublist]

    ntot_2 = len(np.unique(flat_list))


    return ntot_2


# Number of links up to 3 links of distance that connects 2 metal atoms
def get_mmlinks(GG, site, auagposs, upto):
    natoms = len(GG)
    asite = site
    pcont = []
    for i in range(natoms):
        ag_s = i

        sp = nx.shortest_path(G=GG,source=asite,target=ag_s)
        pcont.append(sp)

    pcont = np.array(pcont)

    x = [len(pcont[i]) < upto for i in range(natoms)] # upto 3 links
    x1 = pcont[np.nonzero(x)[0]]

    #Au or Ag nodes from the adsorption site to upto 3 links away
    x2  = [ list(set(auagposs).intersection(set(x1[j])) ) for j in range(len(x1))]

    # flat_list = [item for sublist in x2 for item in sublist]


    count_mml = 0
    for it in x2:
        c = list(itertools.combinations(it, 2)) # All the 2-node combinations
    #     print(c, len(c))
        if len(c)>0:
            for ij in range(0,len(c)):
    #           Finds the distance between 2 nodes in any combination
                splen = nx.shortest_path_length(G=GG, source=c[ij][0], target=c[ij][1])
    #           If the distance = 1, there is a connection between two metal atoms
                if splen == 1:
                    count_mml+=1 # Counts the number of such combinations

    return count_mml

def get_dmat(mol, dict_corad):
    natoms = len(mol)
    dmatrix = np.zeros((natoms,natoms))
    pos = mol[:,1:]
    s = mol[:,0]
    for j in range(natoms):

        d = utils.dist_between(p1 = pos, p2 = pos[j])
#         s = np.array(stf.symbols)

        for i in range(0,natoms):
            dideal = dict_ari[s[j]] + dict_ari[ s[i] ] +.25
            if (d[i] < dideal and d[i]>0):
                dmatrix[j, i ] = 1

    return dmatrix
