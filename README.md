# ugropy

ugropy is a `Python` library to obtain the UNIFAC's subgroups from 
both the name or the SMILES representation of a molecule. If the name is 
given, the library uses the
[PubChemPy](https://github.com/mcs07/PubChemPy) library to obtain the SMILES
representation from PubChem. In both cases, ugropy uses the 
[RDKit](https://github.com/rdkit/rdkit) to search the functional groups in the
molecule.

ugropy is in an early development stage, leaving issues of examples of molecules that ugropy fails solving the UNIFAC's groups is very helpful.

## Example of use
Get UNIFAC groups from the molecule's name:

```python
from ugropy import Groups


hexane = Groups("hexane")
print(hexane.unifac_groups)
```

    {'CH3': 2, 'CH2': 4}

Get UNIFAC groups from molecule's SMILES:

```python
propanol = Groups("CCCO", "smiles")
print(propanol.unifac_groups)
```

    {'CH3': 1, 'CH2': 2, 'OH': 1}

## Installation
At the moment ugropy is not uploaded in PyPI (will be soon).

```
pip install git+https://github.com/SalvadorBrandolin/ugropy.git
```
## Referece

[1] http://www.ddbst.com/published-parameters-unifac.html
