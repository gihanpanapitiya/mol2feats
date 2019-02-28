# mol2feats

A collection of functions to create features for molecular systems.
![alt text](https://github.com/gihanpanapitiya/mol2feats/tree/master/mol2feats/molf.png)



# Installation
(1) Download mol2feats directory </br>
(2) Inside mol2feats directory, type,
```bash
    sudo pip3 install -e .
```

# Usage
## Features based on counts of molecular building blocks
#### Atoms (A) in layers

```python
from mol2feats import counts
from mol2feats import utils

# Read the molecular file
mol = utils.read_mol("benzene.xyz")

# Input:: 
# mol: molecule considered, 
# site: index of the atom around which the boundaries of the layers are centered,
# startl: 1st boundary of the layer goeas through atom identified by the index=startl,
# endll: 2nd boundary of the layer goeas through atom identified by the index=endl,

atoms_in_layers = counts.atoms_inlayers(mol=mol,site=0,startl=0, endl=6)

```

#### Bonds (A-B) in layers
```python

# This function can be used to automatically find bond lenghts
# Input:: 
# Symbols and positions of atoms in the molecule  

blengths = utils.get_bond_lengths(mol_ref = mol)

# Output:: 
A dictionary containing different atom pairs and the corresponding bond lenghts
# blengths = {1.0: {1.0: 2.4788565912533138, 6.0: 1.0830000000000002}, 6.0: {6.0: 1.3959999999999999}}


bonds_in_lyers = counts.bonds_inlayers(mol = mol, site = 0, startl = 0, endl = 6, dB1B2 = blengths)
```


#### A-B-C fragments in layers (A,B and C are different atom types)

```python

# This function calculates the counts of different A-B-C fragments
# Input:: 
# mol: molecule considered, 
# site: index of the atom around which the boundaries of the layers are centered,
# startl: 1st boundary of the layer goeas through atom identified by the index=startl,
# endll: 2nd boundary of the layer goeas through atom identified by the index=endl,
# dB1B2: the dictionary containing the reference bond lengths

abc = counts.abc_in_layers(mol = mol, site = 0, startl=0, endl = 6, dB1B2 = blengths)

# Output:: 
# 1st columns counts of A-B-C
# Columns 2-4: A,B and C
#[[1, 1.0, 1.0, 1.0],
# [4, 1.0, 1.0, 6.0],
# [0, 1.0, 6.0, 1.0],
# [4, 1.0, 6.0, 6.0],
# [4, 6.0, 1.0, 1.0],
# [0, 6.0, 1.0, 6.0],
# [4, 6.0, 6.0, 1.0],
# [1, 6.0, 6.0, 6.0]]

```

## Features based on the distance between atoms

#### Clustering of atoms
```python
# This example calcukates the clustering of all the C atoms with respect to the 
# center of the molecule in Benzene

import numpy as np

pos = mol[mol[:,0]==6][:,1:]
site = np.mean(pos)

clustering, centroid = dist.get_cluster_from_rs(pos=pos, site=site)

```

## Features based on the graphical representation of the molecule


#### Number of 'A' atoms with a path length of 'X' from a specific site
``` Python
from mol2feats import network
import networkx as nx

# Dictionary of covalent radii  
dict_corad ={1:.31, 6:.77}
dmat = network.get_dmat(mol=mol, dict_corad=dict_corad)

# Create the graph
ndmat = np.transpose(np.nonzero(dmat))
G = nx.Graph()
G.add_edges_from(ndmat)
mol = utils.read_mol("benzene2.xyz")

# This code fragment finds the number of S atoms withins a path length of 3 from 
# atomic site specified by 'site' parameter

S_sites = np.where(mol[:,0]==16)[0]
numS = network.get_nlinks_atom(atom_locs = S_sites, site=0, GG=G, upto=5)

```

# Reference
Machine-Learning Prediction of CO Adsorption in Thiolated, Ag-Alloyed Au Nanoclusters</br> 
J. Am. Chem. Soc., 2018, 140 (50), pp 17508–17514 </br>
DOI: 10.1021/jacs.8b08800 </br>
https://pubs.acs.org/doi/10.1021/jacs.8b08800
