"""constants module.

Attributes
----------
R : float
    Gas constant [J/mol/K]
ureg : pint.UnitRegistry
    Unit registry of Pint library
"""

from pathlib import Path

from pint import UnitRegistry


# Unit registry
ureg = UnitRegistry()


# constants.py path
_here = Path(__file__).parent

# CSVs path
_csvs = _here / "groupscsv"

# Gas constant [J/mol/K]
R = 8.31446261815324
