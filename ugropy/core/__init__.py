"""Core module.

FragmentationModel subgroups detection functions.
"""

from .checks import (
    check_can_fit_atoms,
    check_has_composed,
    check_has_composed_overlapping,
    check_has_hiden,
    check_has_molecular_weight_right,
)
from .composed import correct_composed
from .draw_molecule import draw
from .fit_atoms_indexes import fit_atoms
from .get_model_groups import get_groups
from .get_rdkit_object import instantiate_mol_object
from .problematics import correct_problematics


__all__ = [
    "check_can_fit_atoms",
    "check_has_composed",
    "draw",
    "check_has_hiden",
    "check_has_molecular_weight_right",
    "check_has_composed_overlapping",
    "correct_composed",
    "fit_atoms",
    "get_groups",
    "instantiate_mol_object",
    "correct_problematics",
]
