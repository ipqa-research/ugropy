"""Core module.

FragmentationModel subgroups detection functions.
"""

from .checks import check_atoms_fragments_presence

from .get_rdkit_object import instantiate_mol_object


__all__ = [
    check_atoms_fragments_presence,
    "instantiate_mol_object",
]
