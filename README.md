![logo](https://github.com/ipqa-research/ugropy/blob/main/logo.png?raw=true)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ipqa-research/ugropy/blob/main/docs/source/tutorial/easy_way.ipynb)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://tldrlegal.com/license/mit-license)
![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue)
[![Docs](https://img.shields.io/badge/docs%20-%20green?style=flat&label=Sphinx&link=https%3A%2F%2Fipqa-research.github.io%2Fugropy%2Findex.html)](https://salvadorbrandolin.github.io/ugropy/)
[![PyPI
version](https://badge.fury.io/py/ugropy.svg)](https://badge.fury.io/py/ugropy)
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
    

`ugropy` is a `Python` library to obtain subgroups from different thermodynamic
group contribution models using both the name or the SMILES representation of a
molecule. If the name is given, the library uses the
[PubChemPy](https://github.com/mcs07/PubChemPy) library to obtain the SMILES
representation from PubChem. In both cases, `ugropy` uses the
[RDKit](https://github.com/rdkit/rdkit) library to search the functional groups
in the molecule.

`ugropy` is tested for `Python` 3.10, 3.11, 3.12, and 3.13 on Linux, Windows
and Mac OS.

You can access the documentation here: [https://ipqa-research.github.io/ugropy/](https://ipqa-research.github.io/ugropy/)

# Try ugropy now
You can try `ugropy` without installing it by clicking on the Colab badge.

You can install `ugropy` by:

```shell
pip install ugropy
```

# Models implemented

## Gibbs / EoS models
- Classic liquid-vapor UNIFAC
- Predictive Soave-Redlich-Kwong (PSRK)
- Dortmund (modified UNIFAC)

## Property estimators
- Joback
- Abdulelah-Gani (beta)

# Writers
`ugropy` allows you to convert the obtained functional groups or estimated
properties to the input format required by the following thermodynamic
libraries:

- [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl)
- [Thermo](https://github.com/CalebBell/thermo)
- [yaeos (Fortran)](https://github.com/ipqa-research/yaeos)


# Example of use
Here is a little taste of `ugropy`, please, check the full tutorial
[here](https://ipqa-research.github.io/ugropy/tutorial/tutorial.html) to see
all it has to offer!

Get groups from the molecule's name:


```python
from ugropy import Groups


hexane = Groups("hexane")

print(hexane.unifac.subgroups)
print(hexane.psrk.subgroups)
print(hexane.dortmund.subgroups)
print(hexane.joback.subgroups)
print(hexane.agani.primary.subgroups)
```

    {'CH3': 2, 'CH2': 4}
    {'CH3': 2, 'CH2': 4}
    {'CH3': 2, 'CH2': 4}
    {'-CH3': 2, '-CH2-': 4}
    {'CH3': 2, 'CH2': 4}

Get groups from molecule's SMILES:

```python
propanol = Groups("CCCO", "smiles")

print(propanol.unifac.subgroups)
print(propanol.psrk.subgroups)
print(propanol.dortmund.subgroups)
print(propanol.joback.subgroups)
print(propanol.agani.primary.subgroups)
```

    {'CH3': 1, 'CH2': 2, 'OH': 1}
    {'CH3': 1, 'CH2': 2, 'OH': 1}
    {'CH3': 1, 'CH2': 2, 'OH (P)': 1}
    {'-CH3': 1, '-CH2-': 2, '-OH (alcohol)': 1}
    {'CH3': 1, 'CH2': 2, 'OH': 1}

Estimate properties with the Joback and Abdulelah-Gani models!

```python
limonene = Groups("limonene")

print(limonene.joback.subgroups)
print(f"{limonene.joback.critical_temperature} K")
print(f"{limonene.joback.vapor_pressure(176 + 273.15)} bar")
```

    {'-CH3': 2, '=CH2': 1, '=C<': 1, 'ring-CH2-': 3, 'ring>CH-': 1, 'ring=CH-': 1, 'ring=C<': 1}
    657.4486692170663 kelvin
    1.0254019428522743 bar

```python
print(limonene.agani.primary.subgroups)
print(limonene.agani.secondary.subgroups)
print(limonene.agani.tertiary.subgroups)
print(f"{limonene.agani.critical_temperature}")
print(limonene.agani.molecular_weight / limonene.agani.liquid_molar_volume)
```

    {'CH3': 2, 'CH2=C': 1, 'CH2 (cyclic)': 3, 'CH (cyclic)': 1, 'CH=C (cyclic)': 1}
    {'CH3-CHm=CHn (m,n in 0..2)': 1, '(CHn=C)cyc-CH3 (n in 0..2)': 1, 'CHcyc-C=CHn (n in 1..2)': 1}
    {}
    640.1457030826214 kelvin
    834.8700605718585 gram / liter

Visualize your results! (The next code creates the `ugropy` logo)

```Python
mol = Groups("CCCC1=C(COC(C)(C)COC(=O)OCC)C=C(CC2=CC=CC=C2)C=C1", "smiles")

mol.unifac.draw(
    title="ugropy",
    width=800,
    height=450,
    title_font_size=50,
    legend_font_size=14
)
```

`ugropy` can obtain multiple solutions, even nonoptimal ones if desired. For
example:

```python
from ugropy import unifac


solutions = unifac.get_groups(
    "9,10-dihydroanthracene",
    search_multiple_solutions=True,
    search_nonoptimal=True
)

for sol in solutions:
    print(sol.subgroups)
```

```
{'ACH': 8, 'AC': 2, 'ACCH2': 2}
{'CH2': 1, 'ACH': 8, 'AC': 3, 'ACCH2': 1}
{'CH2': 2, 'ACH': 8, 'AC': 4}
```

Write down the [Clapeyron.jl](https://github.com/ClapeyronThermo/Clapeyron.jl)
.csv input files.

```python
from ugropy import writers

names = ["limonene", "adrenaline", "Trinitrotoluene"]

grps = [Groups(n) for n in names]

# Write the csv files into a database directory
writers.to_clapeyron(
    molecules_names=names,
    unifac_groups=[g.unifac.subgroups for g in grps],
    psrk_groups=[g.psrk.subgroups for g in grps],
    joback_objects=[g.joback for g in grps],
    path="database"
)
```
Obtain the [Caleb Bell's Thermo](https://github.com/CalebBell/thermo) subgroups

```python
from ugropy import unifac

names = ["hexane", "ethanol"]

grps = [unifac.get_groups(n) for n in names]

[writers.to_thermo(g.subgroups, unifac) for g in grps]
```

```
[{1: 2, 2: 4}, {1: 1, 2: 1, 14: 1}]
```
