from mol2feats import utils
import numpy as np
# import itertools
# from itertools import combinations
from itertools import product
# from itertools import combinations_with_replacement


# utils = utils.utils()
def abc_in_layers(mol, site, startl, endl, dB1B2):
    # read the molecular file
    mol = mol
    site = site
    startl = startl
    endl = endl
    dB1B2 = dB1B2


    # pos and sym sorted in the ascending order of the distance from site
    pos,sym,_ = utils.get_sorted_st(mol=mol, site=site)

    # nearest neighbor pos and sym
    lp = pos[startl:endl]
    ls = sym[startl:endl]

    typs = np.sort(np.unique(mol[:,0]))
    typ_comb = list(product(typs, repeat=3))
#         bond_typs = list(combinations_with_replacement(typs, 2))

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
        # lp = lp


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
                            d1 = utils.dist_between_two(lp[lpos[k]], lp[cpos[i]])
                            d2 = utils.dist_between_two(lp[rpos[j]], lp[cpos[i]])

                            nnH = np.logical_and(d1 < dB1B2[el[0]][el[1]] +.2, d1 > dB1B2[el[0]][el[1]] -.2 )
                            nnS = np.logical_and(d2 < dB1B2[er[0]][er[1]] +.2, d2 > dB1B2[er[0]][er[1]] -.2 )
                            if np.logical_and(nnH,nnS):
                                count+=1
    #                             print(d1,d2, nnH,nnS)
                    else:
                            d1 = utils.dist_between_two(lp[lpos[k]], lp[cpos[i]])
                            d2 = utils.dist_between_two(lp[rpos[j]], lp[cpos[i]])

                            nnH = np.logical_and(d1 < dB1B2[el[0]][el[1]] +.2, d1 > dB1B2[el[0]][el[1]] -.2 )
                            nnS = np.logical_and(d2 < dB1B2[er[0]][er[1]] +.2, d2 > dB1B2[er[0]][er[1]] -.2 )
                            if np.logical_and(nnH,nnS):
                                count+=1
        count_type.append([count]+ list(typ_comb[ix]) )

    return count_type


# Atom count
def atoms_inspheres(mol, site, startl, endl):

    mol = mol
    site = site
    startl = startl
    endl = endl
    # dB1B2 = dB1B2

    f2 = utils()


    # pos and sym sorted in the ascending order of the distance from site
    pos,sym, _ = f2.get_sorted_st(mol=mol, site=site)

    # nearest neighbor pos and sym
    lp = pos[startl:endl]
    ls = sym[startl:endl]

    typs = np.sort(np.unique(mol[:,0]))

    typ_count = []
    for i in range(len(typs)):
#         print (typs[i], np.count_nonzero(ls ==typs[i]))
        typ_count.append([typs[i], np.count_nonzero(ls ==typs[i])])

    return typ_count

# Bond count
def bonds_inspheres(mol, site, startl, endl):

    mol = mol
    site = site
    startl = startl
    endl = endl
    # dB1B2 = dB1B2

    f2 = utils()


    # pos and sym sorted in the ascending order of the distance from site
    pos,sym, _ = f2.get_sorted_st(mol=mol, site=site)

    # nearest neighbor pos and sym
    lp = pos[startl:endl]
    ls = sym[startl:endl]

    typs = np.sort(np.unique(mol[:,0]))
    typ_comb = list(itertools.combinations(typs,2))
    # typ_comb = list(product(typs, repeat=2))
    for i in typs:
        it = float(i)
        typ_comb.append(tuple([it, it ]))

    atm_pos_a = []
    for atom in typs:
        atm_pos = np.where(ls == atom)[0]
        atm_pos_a.append(atm_pos)

    typ_pos = {}
    for i in range(len(typs)):

        typ_pos[typs[i]] = atm_pos_a[i]

    count_type = []
    for ix in range(len(typ_comb)):

        l = typ_comb[ix][0]
        r = typ_comb[ix][1]

        lpos = typ_pos[l]
        rpos = typ_pos[r]

        el = [l,r]
        el = sorted(el)

        count = 0

        for i in range(len(lpos)):
            for j in range(i,len(rpos)):
    #             print(i,j)

                d1 = f2.dist_between_two(lp[lpos[i]], lp[rpos[j]])
                nnH = np.logical_and(d1 < dB1B2[el[0]][el[1]] +.2, d1 > dB1B2[el[0]][el[1]] -.2 )
                if nnH:
#                     print(nnH)
                    count+=1
        count_type.append([count]+ list(typ_comb[ix]) )

    return count_type
