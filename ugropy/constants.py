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
here = Path(__file__).parent.resolve()

R = 8.31446261815324

# Dataframes
# =============================================================================
# UNIFAC
# =============================================================================
with open(f"{here}/groupscsv/unifac/unifac_subgroups.csv", mode="r") as f:
    unifac_subgroups = pd.read_csv(f, sep="|", index_col="group", comment="?")

with open(f"{here}/groupscsv/unifac/unifac_matrix.csv", mode="r") as f:
    unifac_matrix = pd.read_csv(f, sep="|", index_col="group", comment="?")

with open(f"{here}/groupscsv/unifac/ch2_hideouts.csv", mode="r") as f:
    unifac_ch2_hideouts = pd.read_csv(f, index_col="group", comment="?").index

with open(f"{here}/groupscsv/unifac/ch_hideouts.csv", mode="r") as f:
    unifac_ch_hideouts = pd.read_csv(f, index_col="group", comment="?").index

# Problematics
with open(f"{here}/groupscsv/problematic_structures.csv", mode="r") as f:
    problematic_structures = pd.read_csv(
        f, sep="|", index_col="smarts", comment="?"
    )


# =============================================================================
# PSRK
# =============================================================================
with open(f"{here}/groupscsv/psrk/psrk_subgroups.csv", mode="r") as f:
    psrk_subgroups = pd.read_csv(f, sep="|", index_col="group", comment="?")

with open(f"{here}/groupscsv/psrk/psrk_matrix.csv", mode="r") as f:
    psrk_matrix = pd.read_csv(f, sep="|", index_col="group", comment="?")

with open(f"{here}/groupscsv/psrk/ch2_hideouts.csv", mode="r") as f:
    psrk_ch2_hideouts = pd.read_csv(f, index_col="group", comment="?").index

with open(f"{here}/groupscsv/psrk/ch_hideouts.csv", mode="r") as f:
    psrk_ch_hideouts = pd.read_csv(f, index_col="group", comment="?").index


# =============================================================================
# Dortmund
# =============================================================================
with open(f"{here}/groupscsv/dortmund/dortmund_subgroups.csv", mode="r") as f:
    dortmund_subgroups = pd.read_csv(
        f, sep="|", index_col="group", comment="?"
    )

with open(f"{here}/groupscsv/dortmund/dortmund_matrix.csv", mode="r") as f:
    dortmund_matrix = pd.read_csv(f, sep="|", index_col="group", comment="?")

# with open(f"{here}/groupscsv/dortmund/ch2_hideouts.csv", mode="r") as f:
#     dortmund_ch2_hideouts = pd.read_csv(
#         f, index_col="group", comment="?"
#     ).index

# with open(f"{here}/groupscsv/dortmund/ch_hideouts.csv", mode="r") as f:
#     dortmund_ch_hideouts = pd.read_csv(
# f, index_col="group", comment="?").index


# =============================================================================
# Joback
# =============================================================================
with open(f"{here}/groupscsv/joback/joback_subgroups.csv", mode="r") as f:
    joback_subgroups = pd.read_csv(f, sep="|", index_col="group", comment="?")

with open(f"{here}/groupscsv/joback/joback_matrix.csv", mode="r") as f:
    joback_matrix = pd.read_csv(f, sep="|", index_col="group", comment="?")

with open(f"{here}/groupscsv/joback/joback_ch2_hideouts.csv", mode="r") as f:
    joback_ch2_hideouts = pd.read_csv(f, index_col="group", comment="?").index

with open(f"{here}/groupscsv/joback/joback_ch_hideouts.csv", mode="r") as f:
    joback_ch_hideouts = pd.read_csv(f, index_col="group", comment="?").index

with open(
    f"{here}/groupscsv/joback/joback_problematic_structures.csv", mode="r"
) as f:
    joback_problematics = pd.read_csv(
        f, sep="|", index_col="smarts", comment="?"
    )

with open(
    f"{here}/groupscsv/joback/properties_contributions.csv", mode="r"
) as f:
    joback_properties_contibutions = pd.read_csv(f, sep="|", index_col="group")
