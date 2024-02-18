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


# LV-UNIFAC
_uni = f"{_csvs}/unifac"

unifac = FragmentationModel(
    subgroups=_rd(f"{_uni}/unifac_subgroups.csv", "group"),
    main_groups=_rd(f"{_uni}/unifac_maingroups.csv", "no."),
    problematic_structures=_rd(
        f"{_csvs}/problematic_structures.csv", "smarts"
    ),
    ch2_hideouts=_rd(f"{_uni}/ch2_hideouts.csv", "group").index.to_list(),
    ch_hideouts=_rd(f"{_uni}/ch_hideouts.csv", "group").index.to_list(),
)


# PSRK
_psrk = f"{_csvs}/psrk"

psrk = FragmentationModel(
    subgroups=_rd(f"{_psrk}/psrk_subgroups.csv", "group"),
    main_groups=_rd(f"{_psrk}/psrk_maingroups.csv", "no."),
    problematic_structures=_rd(
        f"{_csvs}/problematic_structures.csv", "smarts"
    ),
    ch2_hideouts=_rd(f"{_psrk}/ch2_hideouts.csv", "group").index.to_list(),
    ch_hideouts=_rd(f"{_psrk}/ch_hideouts.csv", "group").index.to_list(),
)

# Dortmund
_do = f"{_csvs}/dortmund"

dortmund = FragmentationModel(
    subgroups=_rd(f"{_do}/dortmund_subgroups.csv", "group"),
    main_groups=_rd(f"{_do}/dortmund_maingroups.csv", "no."),
    problematic_structures=_rd(f"{_do}/dortmund_problematics.csv", "smarts"),
    ch2_hideouts=_rd(f"{_do}/ch2_hideouts.csv", "group").index.to_list(),
    ch_hideouts=_rd(f"{_do}/ch_hideouts.csv", "group").index.to_list(),
)

# Joback
_jo = f"{_csvs}/joback"

joback = FragmentationModel(
    subgroups=_rd(f"{_jo}/joback_subgroups.csv", "group"),
    problematic_structures=_rd(f"{_jo}/joback_problematics.csv", "smarts"),
)
