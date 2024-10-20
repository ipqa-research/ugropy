"""Classic Liquid Vapor UNIFAC FragmentationModel implementation.

Import and use the UNIFAC FragmentationModel with:

.. code-block:: python

    from ugropy import unifac

    # Get groups from molecule's name
    tol = unifac.get_groups("toluene")

    print(tol.subgroups)

    # Get groups from molecule's SMILES
    eth = unifac.get_groups("CCO", "smiles")

    print(eth.subgroups)

Attributes
----------
unifac : GibbsModel
    Classic LV-UNIFAC FragmentationModel :cite:p:`ddbst, unifac1, unifac2,
    unifac3, unifac4, unifac5, unifac6`
"""

from ugropy.constants import _csvs
from ugropy.core.frag_classes.gibbs_model.gibbs_model import GibbsModel
from ugropy.models.read_csv import _rd


# =============================================================================
# LV-UNIFAC
# =============================================================================
_uni = _csvs / "unifac"

_uni_sg = _rd(_uni / "unifac_subgroups.csv", "group")
_uni_info = _rd(_uni / "unifac_info.csv", "group")

unifac = GibbsModel(subgroups=_uni_sg, subgroups_info=_uni_info)
