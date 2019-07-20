import numpy as np
import matplotlib.pyplot as plt

import itertools
from itertools import combinations
from itertools import product
from itertools import combinations_with_replacement


# class utils:
#
#     def __init__():
#         print("")

def dist_between(p1,p2):
    return np.sqrt(np.sum( (p1-p2)**2, axis=1))

def dist_between_two(p1,p2):
    return np.sqrt(np.sum( (p1-p2)**2))

def read_mol(f_xyz):
    mol = np.loadtxt(f_xyz, skiprows=1)
    return mol

def get_sorted_st(mol,site):

    pos = mol[:,1:]
    syms = mol[:,0]
    xd = dist_between(p1=pos, p2=pos[site]) # dist from site
    srtd = np.argsort(xd) # syms
    syms_sorted = np.array(syms)[srtd]
    pos_sorted = pos[srtd]
    dist_sorted = xd[srtd]

    return (pos_sorted, syms_sorted, dist_sorted)


# reference molecule to find the bond lengths
def get_bond_lengths(mol_ref):

    mol_ref = mol_ref
    # uniqe atom types in the molecule
    typs = np.sort(np.unique(mol_ref[:,0]))

    # all the three atom combinations
    typ_comb = list(product(typs, repeat=3))

    bond_typs = list(combinations_with_replacement(typs, 2))

    #  finding the bond distance between atom types
    atm_pos_a = []
    for atom in typs:
        atm_pos = np.where(mol_ref[:,0] == atom)[0]
        atm_pos_a.append(atm_pos)

    typ_pos = {}
    for i in range(len(typs)):
        typ_pos[typs[i]] = atm_pos_a[i]

    dB1B2 = {}
    for i in range(len(typs)):
        dB1B2[typs[i]] = {}



    for it in range(len(bond_typs)):
    #     print(bond_typs[i])
        at1 = bond_typs[it][0]
        at2 = bond_typs[it][1]

        # print (at1, at2)

        x1 = mol_ref[ typ_pos[at1] ]
        x2 = mol_ref[ typ_pos[at2] ]

        dd_a = []
        for  i in range(len(x1)):
            p1 = x1[i][1:]
            for j in range(len(x2)):
                p2 = x2[j][1:]
                dd = dist_between_two(p1=p1, p2=p2)
                if dd != 0:
                    dd_a.append(dd)
                else:
                    dd_a.append(-100) # if the bond does not exist,
        dB1B2[at1][at2] = np.sort(dd_a)[0]

    return dB1B2



def tetrahedron_volume(a, b, c, d):
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6



def convex_hull_volume_bis(pts):
    ch = ConvexHull(pts)
    simplices = np.column_stack((np.repeat(ch.vertices[0], ch.nsimplex), ch.simplices))

    tets = ch.points[simplices]
    tvols  = tetrahedron_volume(tets[:, 0], tets[:, 1], tets[:, 2], tets[:, 3])
    tvols = list(tvols)

    volsum = np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1], tets[:, 2], tets[:, 3]))
    volsum = float(volsum)
    return volsum, tvols
