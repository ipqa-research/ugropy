"""GC-EoS model.
"""

from ugropy.constants import _csvs
from ugropy.core.frag_classes.base.fragmentation_model import (
    FragmentationModel,
)
from ugropy.models.read_csv import _rd

# =============================================================================
# Joback
# =============================================================================
_gc = _csvs / "gceos"

_gc_sg = _rd(_gc / "gceos_subgroups.csv", "group")

gceos = FragmentationModel(_gc_sg)
