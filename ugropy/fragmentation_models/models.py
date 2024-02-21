"""Fragmentation models implemented.

Attributes
----------
unifac : FragmentationModel
    Classic LV-UNIFAC FragmentationModel
psrk: FragmentationModel
    Predictive Soave-Redlich-Kwong FragmentationModel
dortmund: FragmentationModel
    Dortmund UNIFAC FragmentationModel
joback: FragmentationModel
    Joback FragmentationModel
"""

import pandas as pd

from ugropy.constants import _csvs
from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


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


# =============================================================================
# LV-UNIFAC
# =============================================================================
_uni = f"{_csvs}/unifac"

_uni_sg = _rd(f"{_uni}/unifac_subgroups.csv", "group")
_uni_mg = _rd(f"{_uni}/unifac_maingroups.csv", "no.")
_uni_problems = _rd(f"{_csvs}/problematic_structures.csv", "smarts")
_uni_hide = _rd(f"{_uni}/hideouts.csv", "group")

unifac = FragmentationModel(
    subgroups=_uni_sg,
    main_groups=_uni_mg,
    problematic_structures=_uni_problems,
    hideouts=_uni_hide,
)


# =============================================================================
# PSRK
# =============================================================================
_psrk = f"{_csvs}/psrk"

_psrk_sg = _rd(f"{_psrk}/psrk_subgroups.csv", "group")
_psrk_mg = _rd(f"{_psrk}/psrk_maingroups.csv", "no.")
_psrk_problems = _rd(f"{_csvs}/problematic_structures.csv", "smarts")
_psrk_hide = _rd(f"{_psrk}/hideouts.csv", "group")

psrk = FragmentationModel(
    subgroups=_psrk_sg,
    main_groups=_psrk_mg,
    problematic_structures=_psrk_problems,
    hideouts=_psrk_hide,
)

# =============================================================================
# Dortmund
# =============================================================================
_do = f"{_csvs}/dortmund"

_do_sg = _rd(f"{_do}/dortmund_subgroups.csv", "group")
_do_mg = _rd(f"{_do}/dortmund_maingroups.csv", "no.")
_do_problems = _rd(f"{_do}/dortmund_problematics.csv", "smarts")
_do_hide = _rd(f"{_do}/hideouts.csv", "group")

dortmund = FragmentationModel(
    subgroups=_do_sg,
    main_groups=_do_mg,
    problematic_structures=_do_problems,
    hideouts=_do_hide,
)

# =============================================================================
# Joback
# =============================================================================
_jo = f"{_csvs}/joback"

_do_sg = _rd(f"{_jo}/joback_subgroups.csv", "group")
_do_problems = _rd(f"{_jo}/joback_problematics.csv", "smarts")

joback = FragmentationModel(
    subgroups=_do_sg,
    problematic_structures=_do_problems,
)
