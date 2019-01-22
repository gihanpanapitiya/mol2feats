#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

import itertools
from itertools import combinations
from itertools import product



from itertools import combinations_with_replacement


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

    return (pos_sorted, syms_sorted)




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

        print (at1, at2)

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
        dB1B2[at1][at2] = np.sort(dd_a)[0]

    return dB1B2



mol_ref = read_mol(f_xyz="benz.xyz")

dB1B2 = get_bond_lengths(mol_ref=mol_ref)

mol = read_mol(f_xyz="benz.xyz")


def abc_in_layers(mol, site, startl, endl, dB1B2):
    # read the molecular file
    mol = mol
    site = site
    startl = startl
    endl = endl
    dB1B2 = dB1B2



    # pos and sym sorted in the ascending order of the distance from site
    pos,sym = get_sorted_st(mol=mol, site=site)

    # nearest neighbor pos and sym
    lp = pos[startl:endl]
    ls = sym[startl:endl]



    atm_pos_a = []
    for atom in typs:
        atm_pos = np.where(ls == atom)[0]
        atm_pos_a.append(atm_pos)

    typ_pos = {}
    for i in range(len(typs)):

        typ_pos[typs[i]] = atm_pos_a[i]


#     len(typ_comb)
    count_type = []
    for ix in range(len(typ_comb)):
        l = typ_comb[ix][0]
        c = typ_comb[ix][1]
        r = typ_comb[ix][2]

    #     print(l,c,r)

        lpos = typ_pos[l]
        cpos = typ_pos[c]
        rpos = typ_pos[r]
        lp = lp


        # x = three_ABC(c=c,l=l,r=r,cpos=cpos,lpos=lpos,rpos=rpos,lp=lp)

        el = [l,c]
        el = sorted(el)
        er = [r,c]
        er = sorted(er)

        count = 0


        for i in range(len(cpos)):
            for j in range(len(rpos)):
                for k in range(len(lpos)):
                    if l==r:
                        if k>j:
                            d1 = dist_between_two(lp[lpos[k]], lp[cpos[i]])
                            d2 = dist_between_two(lp[rpos[j]], lp[cpos[i]])

                            nnH = np.logical_and(d1 < dB1B2[el[0]][el[1]] +.2, d1 > dB1B2[el[0]][el[1]] -.2 )
                            nnS = np.logical_and(d2 < dB1B2[er[0]][er[1]] +.2, d2 > dB1B2[er[0]][er[1]] -.2 )
                            if np.logical_and(nnH,nnS):
                                count+=1
    #                             print(d1,d2, nnH,nnS)
                    else:
                            d1 = dist_between_two(lp[lpos[k]], lp[cpos[i]])
                            d2 = dist_between_two(lp[rpos[j]], lp[cpos[i]])

                            nnH = np.logical_and(d1 < dB1B2[el[0]][el[1]] +.2, d1 > dB1B2[el[0]][el[1]] -.2 )
                            nnS = np.logical_and(d2 < dB1B2[er[0]][er[1]] +.2, d2 > dB1B2[er[0]][er[1]] -.2 )
                            if np.logical_and(nnH,nnS):
                                count+=1
        count_type.append([count, typ_comb[ix]])

    return count_type
