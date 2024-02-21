"""constants module.

Attributes
----------
R : float
    Gas constant [J/mol/K]
joback_properties : pandas.DataFrame
    Joback subsgroups' contribution to properties.
"""

from pathlib import Path

import pandas as pd


# constants.py path
_here = Path(__file__).parent.resolve()

# CSVs path
_csvs = f"{_here}/groupscsv"

# Gas constant [J/mol/K]
R = 8.31446261815324

# Joback Properties contributions
with open(f"{_csvs}/joback/properties_contrib.csv", mode="r") as f:
    joback_properties = pd.read_csv(f, sep="|", index_col="group", comment="?")
