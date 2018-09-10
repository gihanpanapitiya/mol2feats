# Number of links up to 3 links of distance that connects 2 metal atoms
def get_mmlinks(GG, adsite, auagposs):
    asite = adsite
    pcont = []
    for i in range(natoms):    
        ag_s = i

        sp = nx.shortest_path(G=GG,source=asite,target=ag_s)
        pcont.append(sp)

    pcont = np.array(pcont)

    x = [len(pcont[i]) < 5 for i in range(natoms)] # upto 3 links
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
