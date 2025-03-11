"""Dortmund mod-UNIFAC FragmentationModel implementation.

Import and use the Dortmund FragmentationModel with:

.. code-block:: python

    from ugropy import dortmund

    # Get groups from molecule's name
    tol = dortmund.get_groups("toluene")

    print(tol.subgroups)

    # Get groups from molecule's SMILES
    eth = dortmund.get_groups("CCO", "smiles")

    print(eth.subgroups)

Attributes
----------
dortmund : GibbsModel
    Dortmund mod-UNIFAC FragmentationModel :cite:p:`dor1, dor2, dor3, dor4,
    dor5, dor6, dor7, dor8, dor9, dor10, dor11`
"""

from ugropy.constants import _csvs
from ugropy.core.frag_classes.gibbs_model.gibbs_model import GibbsModel
from ugropy.models.read_csv import _rd


# =============================================================================
# Dortmund UNIFAC
# =============================================================================
_dor = _csvs / "dortmund"

_dor_sg = _rd(_dor / "dortmund_subgroups.csv", "group")
_dor_info = _rd(_dor / "dortmund_info.csv", "group")

dortmund = GibbsModel(
    subgroups=_dor_sg, subgroups_info=_dor_info, calculate_r_q=False
)
