# mol2feats

A collection of functions to create features for molecular systems.


# Installation
(1) Download mol2feats directory </br>
(2) Inside mol2feats directory, type,
```bash
    sudo pip3 install -e .
```

# Usage
## A-B-C fragments in layers (A,B and C are different atom types)

```python
from mol2feats import counts
from mol2feats import utils

# Read the molecular file
mol = utils.read_mol("benzene.xyz")

# This function can be used to automatically find bond lenghts
# Input: Symbols and positions of atoms in the molecule  
blengths = utils.get_bond_lengths(mol_ref=mol)


# For Benzene 'get_bond_lengths' returns a dictionary containing different 
# atom pairs and the corresponding bond lenghts
blengths = {1.0: {1.0: 2.4788565912533138, 6.0: 1.0830000000000002},
 6.0: {6.0: 1.3959999999999999}}

abc = counts.abc_in_layers(mol=mol, site=1,startl=0,endl=10,dB1B2=blengths)

Output:
[[2, 1.0, 1.0, 1.0],
 [6, 1.0, 1.0, 6.0],
 [0, 1.0, 6.0, 1.0],
 [8, 1.0, 6.0, 6.0],
 [6, 6.0, 1.0, 1.0],
 [0, 6.0, 1.0, 6.0],
 [8, 6.0, 6.0, 1.0],
 [6, 6.0, 6.0, 6.0]]

```
