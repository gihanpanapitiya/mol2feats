# mol2feats

A collection of functions to create features for molecular systems.


# Installation
(1) Download mol2feats directory </br>
(2) Inside mol2feats directory, type,
```bash
    sudo pip3 install -e .
```

# Usage
## Features based on counts of molecular building blocks
### (1) Atoms (A) in layers

```python
from mol2feats import counts
from mol2feats import utils

# Read the molecular file
mol = utils.read_mol("benzene.xyz")

# Input:: mol: molecule considered, 
# site: index of the atom around which the boundaries of the layers are centered,
# startl: 1st boundary of the layer goeas through atom identified by the index=startl,
# endll: 2nd boundary of the layer goeas through atom identified by the index=endl,

atoms_in_layers = counts.atoms_inlayers(mol=mol,site=0,startl=0, endl=6)

```

### (2) Bonds (A-B) in layers
```python

# This function can be used to automatically find bond lenghts
# Input:: Symbols and positions of atoms in the molecule  

blengths = utils.get_bond_lengths(mol_ref = mol)

# Output:: A dictionary containing different atom pairs and the corresponding bond lenghts
# blengths = {1.0: {1.0: 2.4788565912533138, 6.0: 1.0830000000000002}, 6.0: {6.0: 1.3959999999999999}}


bonds_in_lyers = counts.bonds_inlayers(mol = mol, site = 0, startl = 0, endl = 6, dB1B2 = blengths)
```


### (1) A-B-C fragments in layers (A,B and C are different atom types)

```python

# This function calculates the counts of different A-B-C fragments
# Input:: mol: molecule considered, site: index of the atom around which the boundaries of the layers are centered,
# startl: 1st boundary of the layer goeas through atom identified by the index=startl,
# endll: 2nd boundary of the layer goeas through atom identified by the index=endl,
# dB1B2: the dictionary containing the reference bond lengths

abc = counts.abc_in_layers(mol = mol, site = 0, startl=0, endl = 6, dB1B2 = blengths)

# Output:: 1st columns counts of A-B-C,
# columns 2-4: A,B and C
#[[1, 1.0, 1.0, 1.0],
# [4, 1.0, 1.0, 6.0],
# [0, 1.0, 6.0, 1.0],
# [4, 1.0, 6.0, 6.0],
# [4, 6.0, 1.0, 1.0],
# [0, 6.0, 1.0, 6.0],
# [4, 6.0, 6.0, 1.0],
# [1, 6.0, 6.0, 6.0]]

```



