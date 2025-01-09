"""AbdulelahGani Tertiary Structures FragmentationModel implementation.

Import and use the AbdulelahGani Tertiary Structures FragmentationModel with:

.. code-block:: python

    from ugropy import abdulelah_gani_t

    # Get groups from molecule's name
    tol = abdulelah_gani_t.get_groups("toluene")

    print(tol.subgroups)

    # Get groups from molecule's SMILES
    eth = abdulelah_gani_t.get_groups("CCO", "smiles")

    print(eth.subgroups)

Attributes
----------
abdulelah_gani_t: AbdulelahGaniPSTModel
    AbdulelahGaniPSTModel FragmentationModel :cite:p:`gani`
"""

from ugropy.constants import _csvs
from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_pst import (
    AbdulelahGaniPSTModel,
)
from ugropy.models.read_csv import _rd


# =============================================================================
# Abdulelah Gani Tertiary Structures FragmentationModel
# =============================================================================
_ag = _csvs / "abdulelah_gani" / "tertiary"

_ag_sg = _rd(_ag / "tertiary.csv", "group")
_ag_info = _rd(_ag / "info.csv", "group")

abdulelah_gani_t = AbdulelahGaniPSTModel(
    subgroups=_ag_sg,
    subgroups_info=_ag_info,
    allow_overlapping=True,
    allow_free_atoms=True,
)
