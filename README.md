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

mol = utils.read_mol("../benzene.xyz")
blengths = utils.get_bond_lengths(mol_ref=mol)

counts.abc_in_layers(mol=mol, site=1,startl=0,endl=10,dB1B2=blengths)

```
