"""Core module.

FragmentationModel subgroups detection functions.
"""

from .checks import (
    check_has_overlapping_groups
)

from .get_rdkit_object import instantiate_mol_object


__all__ = [
    "check_has_overlapping_groups",
    "instantiate_mol_object",
]
