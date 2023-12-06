<p align="center">
  <img src="logo.png" alt="logo" width="300" height="350">
</p>

ugropy is a `Python` library to obtain the UNIFAC's subgroups from both the
name or the SMILES representation of a molecule. If the name is given, the
library uses the [PubChemPy](https://github.com/mcs07/PubChemPy) library to
obtain the SMILES representation from PubChem. In both cases, ugropy uses the
[RDKit](https://github.com/rdkit/rdkit) to search the functional groups in the
molecule.

ugropy is in an early development stage, leaving issues of examples of
molecules that ugropy fails solving the UNIFAC's groups is very helpful.

## Models supported v1.0.0
- Classic liquid-vapor UNIFAC
- Predictive Soave-Redlich-Kwong (PSRK)
- Joback

## Example of use
Get UNIFAC groups from the molecule's name:

```python
from ugropy import Groups


hexane = Groups("hexane")

print(hexane.unifac_groups)
print(hexane.psrk_groups)
print(hexane.joback.groups)
```

    {'CH3': 2, 'CH2': 4}
    {'CH3': 2, 'CH2': 4}
    {'-CH3': 2, '-CH2-': 4}

Get UNIFAC groups from molecule's SMILES:

```python
propanol = Groups("CCCO", "smiles")

print(propanol.unifac_groups)
print(propanol.psrk_groups)
print(propanol.joback.groups)
```

    {'CH3': 1, 'CH2': 2, 'OH': 1}
    {'CH3': 1, 'CH2': 2, 'OH': 1}
    {'-CH3': 1, '-CH2-': 2, '-OH (alcohol)': 1}

Estimate properties with the Joback model!

```python
limonene = Groups("limonene")

print(limonene.joback.groups)
print(f"{limonene.joback.critical_temperature} K")
print(f"{limonene.joback.vapor_pressure(176 + 273.15)} bar")
```

    {'-CH3': 2, '=CH2': 1, '=C<': 1, 'ring-CH2-': 3, 'ring>CH-': 1, 'ring=CH-': 1, 'ring=C<': 1}
    657.4486692170663 K
    1.0254019428522743 bar

## Installation
At the moment ugropy is not uploaded in PyPI (will be soon).

```
pip install git+https://github.com/SalvadorBrandolin/ugropy.git
```

## Refereces

[1] http://www.ddbst.com/published-parameters-unifac.html

[2] Joback, K. G., & Reid, R. C. (1987). ESTIMATION OF PURE-COMPONENT
PROPERTIES FROM GROUP-CONTRIBUTIONS. Chemical Engineering Communications,
57(1–6), 233–243. https://doi.org/10.1080/00986448708960487

[3] Joback, K. G. (1989). Designing molecules possessing desired physical
property values [Thesis (Ph. D.), Massachusetts Institute of Technology].
https://dspace.mit.edu/handle/1721.1/14191
