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


def _rd(file_path: str, index_col: str = None) -> pd.DataFrame:
    """Read the models' csv.

    Parameters
    ----------
    file_path : str
        Path to csv file.
    index_col : str, optional
        Name of the index column, by default None.

    Returns
    -------
    pd.DataFrame
        Readed csv.
    """
    with open(file_path, mode="r") as f:
        return pd.read_csv(f, sep="|", index_col=index_col, comment="?")


# constants.py path
_here = Path(__file__).parent.resolve()


# CSVs path
_csvs = f"{_here}/groupscsv"


# Gas constant [J/mol/K]
R = 8.31446261815324


# Gibss Excess models dataframes
# =============================================================================
# UNIFAC
# =============================================================================
unifac_subgroups = _rd(f"{_csvs}/unifac/unifac_subgroups.csv", "group")
unifac_maingroups = _rd(f"{_csvs}/unifac/unifac_maingroups.csv", "no.")
unifac_matrix = _rd(f"{_csvs}/unifac/unifac_matrix.csv", "group")
unifac_ch2_hide = _rd(f"{_csvs}/unifac/ch2_hideouts.csv", "group").index
unifac_ch_hide = _rd(f"{_csvs}/unifac/ch_hideouts.csv", "group").index
unifac_problem = _rd(f"{_csvs}/problematic_structures.csv", "smarts")


# =============================================================================
# PSRK
# =============================================================================
psrk_subgroups = _rd(f"{_csvs}/psrk/psrk_subgroups.csv", "group")
psrk_maingroups = _rd(f"{_csvs}/psrk/psrk_maingroups.csv", "no.")
psrk_matrix = _rd(f"{_csvs}/psrk/psrk_matrix.csv", "group")
psrk_ch2_hide = _rd(f"{_csvs}/psrk/ch2_hideouts.csv", "group").index
psrk_ch_hide = _rd(f"{_csvs}/psrk/ch_hideouts.csv", "group").index
psrk_problem = _rd(f"{_csvs}/problematic_structures.csv", "smarts")


# =============================================================================
# Dortmund
# =============================================================================
dort_subgroups = _rd(f"{_csvs}/dortmund/dortmund_subgroups.csv", "group")
dort_maingroups = _rd(f"{_csvs}/dortmund/dortmund_maingroups.csv", "no.")
dort_matrix = _rd(f"{_csvs}/dortmund/dortmund_matrix.csv", "group")
dort_ch2_hide = _rd(f"{_csvs}/dortmund/ch2_hideouts.csv", "group").index
dort_ch_hide = _rd(f"{_csvs}/dortmund/ch_hideouts.csv", "group").index
dort_problem = _rd(f"{_csvs}/dortmund/dortmund_problematics.csv", "smarts")


# Properties estimators models dataframes
# =============================================================================
# Joback
# =============================================================================
joback_subgroups = _rd(f"{_csvs}/joback/joback_subgroups.csv", "group")
joback_matrix = _rd(f"{_csvs}/joback/joback_matrix.csv", "group")
joback_ch2_hide = _rd(f"{_csvs}/joback/ch2_hideouts.csv", "group").index
joback_ch_hide = _rd(f"{_csvs}/joback/ch_hideouts.csv", "group").index
joback_problem = _rd(f"{_csvs}/joback/joback_problematics.csv", "smarts")
joback_properties = _rd(f"{_csvs}/joback/properties_contrib.csv", "group")
