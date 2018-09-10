dict1 = {47.:'S',79.:'G'}
ref4bonds = ['GGGG', 'GGGS', 'GGSG', 'GSGG', 'SGGG', 'SGGS']

def fouratombonds(p_sorted, s_sorted, l_sorted, strt_n, end_n):
#     n_nn = n_neighs
    ss = s_sorted[strt_n : end_n]
    sp = np.where(np.logical_or(ss==79, ss==47))

    s_sel = ss[sp]
    p_sel = p_sorted[strt_n : end_n][sp]

    l_sel = l_sorted[strt_n : end_n][sp]

    # Au-Au-Au-Au
    bonds = []
#     print(len(p_sel))
    for i in range(0,len(p_sel) ):
    #     if s_sel[i] == 79:
        p1 = p_sel[i]
        s1 = s_sel[i]
        l1 = l_sel[i]
#         print(i, p_sel)
        d = dist_between(p1=p_sel,p2=p1)
        asrt  = np.argsort(d)
        if d[asrt][1] < 3.2:
            p2 = p_sel[asrt][1]
            s2 = s_sel[asrt][1]
            l2 = l_sel[asrt][0:len(p_sel)][1]
            p_sel2 = p_sel[asrt][1:len(p_sel)]
            s_sel2 = s_sel[asrt][1:len(p_sel)]
            l_sel2 = l_sel[asrt][1:len(p_sel)]
            d2 = dist_between(p1=p_sel2,p2=p2)
            asrt2  = np.argsort(d2)
            if d2[asrt2][1] < 3.2:
                p3 = p_sel2[asrt2][1]
                s3 = s_sel2[asrt2][1]
                l3 = l_sel2[asrt2][1]
                p_sel3 = p_sel2[asrt2][1:len(p_sel2)]
                s_sel3 = s_sel2[asrt2][1:len(p_sel2)]
                l_sel3 = l_sel2[asrt2][1:len(p_sel2)]
                d3 = dist_between(p1=p_sel3,p2=p3)
                asrt3  = np.argsort(d3)
                if d3[asrt3][1] < 3.2:
                    p4 = p_sel3[asrt3][1]
                    s4 = s_sel3[asrt3][1]
                    l4 = l_sel3[asrt3][1]
    #                 p_sel3 = p_sel2[asrt2][1:len(p_sel2)]
    #                 d4 = dist_between(a=p_sel4,p=p4)
    #                 asrt4  = np.argsort(d4)
#         print(s1,s2,s3,s4, d[asrt][1],d2[asrt2][1],d3[asrt3][1], l1,l2,l3,l4 )
        bonds.append(dict1[s1]+dict1[s2]+dict1[s3]+dict1[s4] )
    
    bonds = np.array(bonds)
    u_bonds = []
    for j in range(0,len(ref4bonds) ):
        un = len(np.where(bonds == ref4bonds[j])[0])
        u_bonds.append(un)

    return u_bonds
