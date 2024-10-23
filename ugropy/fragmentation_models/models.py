"""Fragmentation models implemented.

All models can be imported directly with:

.. code-block:: python

   from ugropy import joback, psrk, unifac

You can check the group list and their SMARTS representation with:

.. code-block:: python

   joback.subgroups psrk.subgroups unifac.subgroups

In the case of a PropertiesEstimator like joback, you can check the
contribution of each group to the properties with:

.. code-block:: python

    joback.properties_contributions

Attributes
----------
unifac : GibbsModel
    Classic LV-UNIFAC FragmentationModel :cite:p:`ddbst, unifac1, unifac2,
    unifac3, unifac4, unifac5, unifac6`
psrk: GibbsModel
    Predictive Soave-Redlich-Kwong FragmentationModel :cite:p:`ddbst, psrk1,
    psrk2`
joback: PropertiesEstimator
    Joback FragmentationModel :cite:p:`joback1, joback2`
"""

import pandas as pd

from ugropy.constants import _csvs
from ugropy.fragmentation_models.gibbs_model import GibbsModel
from ugropy.fragmentation_models.prop_estimator import PropertiesEstimator


def _rd(file_path: str, index_col: str = None) -> pd.DataFrame:
    """Read the models' csv.

    Parameters
    ----------
    file_path : str or pathlib
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
_uni = _csvs / "unifac"

_uni_sg = _rd(_uni / "unifac_subgroups.csv", "group")
_uni_mg = _rd(_uni / "unifac_maingroups.csv", "no.")
_uni_info = _rd(_uni / "unifac_info.csv", "group")
_uni_problems = _rd(_csvs / "problematic_structures.csv", "smarts")
_uni_hide = _rd(_uni / "hideouts.csv", "group")

unifac = GibbsModel(
    subgroups=_uni_sg,
    split_detection_smarts=["C5H4N", "C5H3N", "C4H3S", "C4H2S"],
    problematic_structures=_uni_problems,
    hideouts=_uni_hide,
    subgroups_info=_uni_info,
    main_groups=_uni_mg,
)


# =============================================================================
# PSRK
# =============================================================================
_psrk = _csvs / "psrk"

_psrk_sg = _rd(_psrk / "psrk_subgroups.csv", "group")
_psrk_mg = _rd(_psrk / "psrk_maingroups.csv", "no.")
_psrk_info = _rd(_psrk / "psrk_info.csv", "group")
_psrk_problems = _rd(_csvs / "problematic_structures.csv", "smarts")
_psrk_hide = _rd(_psrk / "hideouts.csv", "group")

psrk = GibbsModel(
    subgroups=_psrk_sg,
    split_detection_smarts=["C5H4N", "C5H3N", "C4H3S", "C4H2S"],
    problematic_structures=_psrk_problems,
    hideouts=_psrk_hide,
    subgroups_info=_psrk_info,
    main_groups=_psrk_mg,
)

# =============================================================================
# Dortmund
# =============================================================================
# _do = f"{_csvs}/dortmund"

# _do_sg = _rd(f"{_do}/dortmund_subgroups.csv", "group") _do_mg =
# _rd(f"{_do}/dortmund_maingroups.csv", "no.") _do_problems =
# _rd(f"{_do}/dortmund_problematics.csv", "smarts") _do_hide =
# _rd(f"{_do}/hideouts.csv", "group")
#
# dortmund = FragmentationModel( subgroups=_do_sg, main_groups=_do_mg,
#    problematic_structures=_do_problems, hideouts=_do_hide, )

# =============================================================================
# Joback
# =============================================================================
_jo = _csvs / "joback"

_jo_sg = _rd(_jo / "joback_subgroups.csv", "group")
_jo_problems = _rd(_jo / "joback_problematics.csv", "smarts")
_jo_props = _rd(_jo / "properties_contrib.csv", "group")

joback = PropertiesEstimator(
    subgroups=_jo_sg,
    problematic_structures=_jo_problems,
    properties_contributions=_jo_props,
)


# =============================================================================
# Constantinou and Gani
# =============================================================================
# Primary structures
_cg_p = _csvs / "constantinou_gani" / "primary"

_cg_sg = _rd(_cg_p / "c_g_prymary_subgroups.csv", "group")
_cg_problems = _rd(_cg_p / "cg_problematics.csv", "smarts")
_cg_hide = _rd(_cg_p / "hideouts.csv", "group")
_cg_props = _rd(_cg_p / "properties_prymary_contrib.csv", "group")

constantinou_gani_primary = PropertiesEstimator(
    subgroups=_cg_sg,
    split_detection_smarts=["C5H4N", "C5H3N", "C4H3S", "C4H2S"],
    problematic_structures=_cg_problems,
    hideouts=_cg_hide,
    properties_contributions=_cg_props,
)
