def get_spathag(GG, adsite, auagposs, agposs ):
    zarr = []
    for ix in auagposs:
        z = nx.shortest_path_length(G=GG, source=adsite, target=ix)
#         if z <3:
        if ix in agposs:
            z = 200+z
#             print (z)
        zarr.append(z)
    return zarr
