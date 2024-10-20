"""check module.

The module contains the necessary checks to corroborate the success of the
algorithm to obtain the molecule's FragmentationModel subgroups.

"""

import numpy as np

from rdkit import Chem


def check_atoms_fragments_presence(
    molecule: Chem.rdchem.Mol, fragments: dict
) -> tuple[bool, np.ndarray]:
    """Find overlapped atoms and free atoms.

    Check the detected fragments to find the atoms that appears in more
    than one fragment (overlapping), and the atoms that are not present in
    any fragment (free atoms). Returning two np.ndarray with the indexes of
    the overlapping and free atoms.

    Example of a `fragments` dictionary that not presents overlapping
    atoms. For example N-hexane:

    .. code-block:: python

        {
            'CH3_0': (0,),
            'CH3_1': (5,),
            'CH2_0': (1,),
            'CH2_1': (2,),
            'CH2_2': (3,),
            'CH2_3': (4,)
        }

    Example of a `fragments` dictionary that presents overlapping atoms. For
    exmple toluene:

    .. code-block:: python

        {
            'CH3_0': (0,),
            'ACH_0': (2,),
            'ACH_1': (3,),
            'ACH_2': (4,),
            'ACH_3': (5,),
            'ACH_4': (6,),
            'AC_0': (1,),
            'ACCH3_0': (1, 0)
        }

    Parameters
    ----------
    molecule : Chem.rdchem.Mol
        RDKit molecule object.
    fragments : dict
        Dictionary containing the fragments detected in the molecule. The
        keys are the group names and the values are the indexes of the
        atoms in the group.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Overlapping atoms indexes and free atoms indexes.
    """
    n_atoms = molecule.GetNumAtoms()

    # Count the number of times an atom is in a group. Also find the atoms
    # that are not present in any fragment.
    atoms = np.zeros(n_atoms, dtype=int)

    for indexes in fragments.values():
        np.add.at(atoms, np.array(indexes).flatten(), 1)

    overlapped_atoms = np.argwhere(atoms > 1).flatten()
    free_atoms = np.argwhere(atoms == 0).flatten()

    return overlapped_atoms, free_atoms


# def check_has_molecular_weight_right(
#     self, mol_object: Chem.rdchem.Mol, fragments: dict
# ) -> bool:

#     # rdkit molecular weight
#     rdkit_mw = Descriptors.MolWt(mol_object)

#     # Molecular weight from functional groups
#     mws = self.mol_subgroups.loc[
#         list(fragments.keys()), "molecular_weight"
#     ].to_numpy()

#     func_group_mw = np.dot(mws, list(mol_subgroups.values()))

#     return np.allclose(rdkit_mw, func_group_mw, atol=0.5)
