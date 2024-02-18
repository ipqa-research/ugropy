"""constants module.

Attributes
----------
R : float
    Gas constant [J/mol/K]
unifac_subgroups : pandas.DataFrame
    Classic LV-UNIFAC subgroups with it's SMARTS representation, contribution
    and composed classification.
unifac_matrix : pandas.DataFrame
    Classic LV-UNIFAC contribution matrix.
unifac_ch2_hideouts : pandas.DataFrame
    Classic LV-UNIFAC CH2 hideouts.
unifac_ch_hideouts : pandas.DataFrame
    Classic LV-UNIFAC CH hideouts.
problematic_structures : pandas.DataFrame
    Problematic structures.
psrk_subgroups : pandas.DataFrame
    Classic PSRK subgroups with it's SMARTS representation, contribution and
    composed classification.
psrk_matrix : pandas.DataFrame
    Classic PSRK contribution matrix.
psrk_ch2_hideouts : pandas.DataFrame
    Classic PSRK CH2 hideouts.
psrk_ch_hideouts : pandas.DataFrame
    Classic PSRK CH hideouts.
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
