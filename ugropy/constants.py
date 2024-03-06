"""constants module.

Attributes
----------
R : float
    Gas constant [J/mol/K]
"""

from pathlib import Path


# constants.py path
_here = Path(__file__).parent.resolve()

# CSVs path
_csvs = f"{_here}/groupscsv"

# Gas constant [J/mol/K]
R = 8.31446261815324
