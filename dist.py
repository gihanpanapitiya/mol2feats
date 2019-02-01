from mol2feats import utils
import numpy as np


# utils = utils.utils()
def get_cluster_from_rs(pos, site):
    # f2 = utils()

    pos = pos
    site = site


    if len(pos) == 0:
        cluster_f, centroid = 0,0

    elif len(pos) == 1:
        centroid = np.mean(pos, axis=0)
        cluster_f  = utils.dist_between_two(pos , centroid)

    else:
        centroid = np.mean(pos, axis=0)
        cluster_f  = np.mean(np.sqrt(np.sum((pos - centroid)**2, axis=1)))

    return (cluster_f, centroid)
