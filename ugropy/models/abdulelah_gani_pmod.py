"""AbdulelahGani Primary Structures FragmentationModel implementation.

Import and use the AbdulelahGani Primary Structures FragmentationModel with:

.. code-block:: python

    from ugropy import abdulelah_gani_p

    # Get groups from molecule's name
    tol = abdulelah_gani_p.get_groups("toluene")

    print(tol.subgroups)

    # Get groups from molecule's SMILES
    eth = abdulelah_gani_p.get_groups("CCO", "smiles")

    print(eth.subgroups)

Attributes
----------
abdulelah_gani_p: AbdulelahGaniPrimaryModel
    AbdulelahGaniPrimaryModel FragmentationModel :cite:p:`gani`
"""

from ugropy.constants import _csvs
from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_p import (
    AbdulelahGaniPrimaryModel,
)
from ugropy.models.read_csv import _rd


# =============================================================================
# Abdulelah Gani Primary Structures FragmentationModel
# =============================================================================
_ag = _csvs / "abdulelah_gani" / "primary"

_ag_sg = _rd(_ag / "primary.csv", "group")
_ag_info = _rd(_ag / "info.csv", "group")

abdulelah_gani_p = AbdulelahGaniPrimaryModel(_ag_sg, _ag_info)
