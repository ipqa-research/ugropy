"""AbdulelahGani Secondary Structures FragmentationModel implementation.

Import and use the AbdulelahGani Secondary Structures FragmentationModel with:

.. code-block:: python

    from ugropy import abdulelah_gani_s

    # Get groups from molecule's name
    tol = abdulelah_gani_s.get_groups("toluene")

    print(tol.subgroups)

    # Get groups from molecule's SMILES
    eth = abdulelah_gani_s.get_groups("CCO", "smiles")

    print(eth.subgroups)

Attributes
----------
abdulelah_gani_s: AbdulelahGaniPSTModel
    AbdulelahGaniPSTModel FragmentationModel :cite:p:`gani`
"""

from ugropy.constants import _csvs
from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_pst import (
    AbdulelahGaniPSTModel,
)
from ugropy.models.read_csv import _rd


# =============================================================================
# Abdulelah Gani Secondary Structures FragmentationModel
# =============================================================================
_ag = _csvs / "abdulelah_gani" / "secondary"

_ag_sg = _rd(_ag / "secondary.csv", "group")
_ag_info = _rd(_ag / "info.csv", "group")

abdulelah_gani_s = AbdulelahGaniPSTModel(
    subgroups=_ag_sg,
    subgroups_info=_ag_info,
    allow_overlapping=True,
    allow_free_atoms=True,
)
