"""Core module.

UNIFAC's subgroups detection functions.
"""

from .checks import (
    check_has_composed,
    check_has_hidden_ch2_ch,
    check_has_molecular_weight_right,
    check_has_hidden
)
from .composed import correct_composed
from .get_model_groups import get_groups
from .get_rdkit_object import instantiate_mol_object
from .problematics import correct_problematics


__all__ = [
    "check_has_composed",
    "check_has_hidden_ch2_ch",
    "check_has_molecular_weight_right",
    "correct_composed",
    "get_groups",
    "instantiate_mol_object",
    "correct_problematics",
]
