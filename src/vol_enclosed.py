def get_volagrs(atom_sites, rs_pos, stlcl, n_atoms):

    ags2 = np.append(arr = atom_sites, values=rs_pos)
    agpos = stlcl.positions[ags2]

    if n_atoms == 0:
        volagrs = 0
    elif n_atoms == 1:
        p1, p2 = stlcl.positions[atom_sites], stlcl.positions[rs_pos]
        volagrs = dist_between_two(p1,p2)
    elif n_atoms == 2:
        p1, p2, p3 = stlcl.positions[atom_sites[0]], stlcl.positions[atom_sites[1]], stlcl.positions[rs_pos]
        volagrs = trianglearea([p1,p2,p3])
    elif (n_atoms == 3) and (rs_pos in atom_sites):
        print(atom_sites, rs_pos)
        loc = np.where(rs_pos != atom_sites)
        p1, p2, p3 = stlcl.positions[atom_sites[loc][0]], stlcl.positions[atom_sites[loc][1]],\
        stlcl.positions[rs_pos]
        volagrs = trianglearea([p1,p2,p3])
    else:
        volagrs, _ = convex_hull_volume_bis(agpos)

    return volagrs
