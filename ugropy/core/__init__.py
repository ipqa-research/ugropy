"""Core module.

The main algorithm functions and classes are included in this module.

- ilp_solvers:
    Module with the structure of the integer linear programming (ILP) solvers.
- frag_classes:
    Module with the classes that define the fragmentation models.
- check_atoms_fragments_presence:
    Check function to corroborate if ILP algorithm is necessary to solve a
    molecule.
- instantiate_mol_object:
    Function to instantiate a RDKit molecule object from a SMILES string, name
    or rdkit mol object.
"""

from . import frag_classes, ilp_solvers
from .checks import check_atoms_fragments_presence
from .get_rdkit_object import instantiate_mol_object


__all__ = [
    "ilp_solvers",
    "frag_classes",
    "check_atoms_fragments_presence",
    "instantiate_mol_object",
]
