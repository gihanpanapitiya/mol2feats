def get_cluster(agpos, st_loc):
    st_loc = st_loc
    agpos = agpos
    li = list(set(core).intersection(set(   agpos  )))

    centroid = np.mean(st_loc.positions[li], axis=0)

    all_ag_pos = st_loc.positions[  agpos  ]
    core_ag_pos = st_loc.positions[  li  ]
    

    cluster_f  = np.mean(np.sqrt(np.sum((core_ag_pos - centroid)**2, axis=1)))
    cluster_fall  = np.mean(np.sqrt(np.sum((all_ag_pos - centroid)**2, axis=1)))
    
    return (cluster_f, centroid, cluster_fall)
