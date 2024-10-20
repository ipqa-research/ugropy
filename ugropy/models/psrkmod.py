"""Predictive Soave-Redlich-Kwong (PSRK) FragmentationModel implementation.

Import and use the PSRK FragmentationModel with:

.. code-block:: python

    from ugropy import psrk

    # Get groups from molecule's name
    tol = psrk.get_groups("toluene")

    print(tol.subgroups)

    # Get groups from molecule's SMILES
    eth = psrk.get_groups("CCO", "smiles")

    print(eth.subgroups)

Attributes
----------
psrk: GibbsModel
    Predictive Soave-Redlich-Kwong FragmentationModel :cite:p:`ddbst, psrk1,
    psrk2`
"""

from ugropy.constants import _csvs
from ugropy.core.frag_classes.gibbs_model.gibbs_model import GibbsModel
from ugropy.models.read_csv import _rd


# =============================================================================
# PSRK
# =============================================================================
_psrk = _csvs / "psrk"

_psrk_sg = _rd(_psrk / "psrk_subgroups.csv", "group")
_psrk_info = _rd(_psrk / "psrk_info.csv", "group")

psrk = GibbsModel(subgroups=_psrk_sg, subgroups_info=_psrk_info)
