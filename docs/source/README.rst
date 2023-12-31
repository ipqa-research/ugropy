|Logo|

|Binder| |License| |Python 3.10+| |Documentation Status| |pypi|

ugropy
======

ugropy is a ``Python`` library to obtain the UNIFAC’s subgroups from
both the name or the SMILES representation of a molecule. If the name is
given, the library uses the
`PubChemPy <https://github.com/mcs07/PubChemPy>`__ library to obtain the
SMILES representation from PubChem. In both cases, ugropy uses the
`RDKit <https://github.com/rdkit/rdkit>`__ library to search the
functional groups in the molecule.

ugropy is in an early development stage, leaving issues of examples of
molecules that ugropy fails solving the UNIFAC’s groups is very helpful.

Models supported v1.0.0
=======================

-  Classic liquid-vapor UNIFAC
-  Predictive Soave-Redlich-Kwong (PSRK)
-  Joback

Example of use
==============

You can check the full tutorial
`here <https://ugropy.readthedocs.io/en/latest/tutorial/tutorial.html>`__.

Get UNIFAC groups from the molecule’s name:

.. code:: python

   from ugropy import Groups


   hexane = Groups("hexane")

   print(hexane.unifac_groups)
   print(hexane.psrk_groups)
   print(hexane.joback.groups)

::

   {'CH3': 2, 'CH2': 4}
   {'CH3': 2, 'CH2': 4}
   {'-CH3': 2, '-CH2-': 4}

Get UNIFAC groups from molecule’s SMILES:

.. code:: python

   propanol = Groups("CCCO", "smiles")

   print(propanol.unifac_groups)
   print(propanol.psrk_groups)
   print(propanol.joback.groups)

::

   {'CH3': 1, 'CH2': 2, 'OH': 1}
   {'CH3': 1, 'CH2': 2, 'OH': 1}
   {'-CH3': 1, '-CH2-': 2, '-OH (alcohol)': 1}

Estimate properties with the Joback model!

.. code:: python

   limonene = Groups("limonene")

   print(limonene.joback.groups)
   print(f"{limonene.joback.critical_temperature} K")
   print(f"{limonene.joback.vapor_pressure(176 + 273.15)} bar")

::

   {'-CH3': 2, '=CH2': 1, '=C<': 1, 'ring-CH2-': 3, 'ring>CH-': 1, 'ring=CH-': 1, 'ring=C<': 1}
   657.4486692170663 K
   1.0254019428522743 bar

Write down the Clapeyron.jl .csv input files.

.. code:: python

   from ugropy import writers

   names = ["limonene", "adrenaline", "Trinitrotoluene"]

   grps = [Groups(n) for n in names]

   # Write the csv files into a database directory
   writers.to_clapeyron(
      molecules_names=names,
      unifac_groups=[g.unifac_groups for g in grps],
      psrk_groups=[g.psrk_groups for g in grps],
      joback_objects=[g.joback for g in grps],
      path="./database"
   )

Installation
============

At the moment ugropy is not uploaded in PyPI (will be soon).

::

   pip install git+https://github.com/ipqa-research/ugropy.git

Refereces
=========

[1] http://www.ddbst.com/published-parameters-unifac.html

[2] Joback, K. G., & Reid, R. C. (1987). ESTIMATION OF PURE-COMPONENT
PROPERTIES FROM GROUP-CONTRIBUTIONS. Chemical Engineering
Communications, 57(1–6), 233–243.
https://doi.org/10.1080/00986448708960487

[3] Joback, K. G. (1989). Designing molecules possessing desired
physical property values [Thesis (Ph. D.), Massachusetts Institute of
Technology]. https://dspace.mit.edu/handle/1721.1/14191

.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/ipqa-research/ugropy/main
.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://tldrlegal.com/license/mit-license
.. |Python 3.10+| image:: https://img.shields.io/badge/Python-3.10%2B-blue
.. |Documentation Status| image:: https://readthedocs.org/projects/ugropy/badge/?version=latest
   :target: https://ugropy.readthedocs.io/en/latest/?badge=latest
.. |Logo| image:: https://github.com/ipqa-research/ugropy/blob/main/logo.png?raw=true
   :target: https://github.com/ipqa-research/ugropy
.. |pypi| image:: https://badge.fury.io/py/ugropy.svg
   :target: https://badge.fury.io/py/ugropy
