"""Core module.

UNIFAC's subgroups detection functions.
"""
from .get_groups import get_groups
from .model_getters import (
    get_joback_groups,
    get_psrk_groups,
    get_unifac_groups,
    instantiate_chem_object,
)


__all__ = [
    "get_groups",
    "get_joback_groups",
    "get_psrk_groups",
    "get_unifac_groups",
    "instantiate_chem_object",
]
