def get_cluster_from_rs(agpos, st_loc, rs):
# Agsites = np.where(np.array(st.symbols)=='Ag')[0]
# st_loc = st
# agpos = Agsites
# rs = rs

    st_loc = st_loc
    agpos = agpos
    rs = rs

    li = list(set(core).intersection(set(   agpos  )))

    all_ag_pos = st_loc.positions[  agpos  ] - st_loc.positions[rs]
    core_ag_pos = st_loc.positions[  li  ] - st_loc.positions[rs]



    if len(all_ag_pos) == 0:
        cluster_f,cluster_fall, cluster_fall2,centroid, centroid_all = 0,0,0,0,0

    elif len(all_ag_pos) == 1:
        centroid = np.mean(core_ag_pos, axis=0)
        centroid_all = np.mean(all_ag_pos, axis=0)

        cluster_f  = dist_between_two(core_ag_pos , centroid)
        cluster_fall  = dist_between_two(all_ag_pos , centroid)
#         cluster_fall  = dist_between_two(all_ag_pos , centroid_all) # changing the definition just for 1Ag 
        cluster_fall2  = dist_between_two(all_ag_pos , centroid_all)


    else:
        centroid = np.mean(core_ag_pos, axis=0)
        centroid_all = np.mean(all_ag_pos, axis=0)

        cluster_f  = np.mean(np.sqrt(np.sum((core_ag_pos - centroid)**2, axis=1)))
        cluster_fall  = np.mean(np.sqrt(np.sum((all_ag_pos - centroid)**2, axis=1)))
        cluster_fall2  = np.mean(np.sqrt(np.sum((all_ag_pos - centroid_all)**2, axis=1)))

    return (cluster_f, centroid, cluster_fall,cluster_fall2, centroid_all)
