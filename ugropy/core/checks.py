"""check module.

The module contains the necessary checks to corroborate the success of the
algorithm to obtain the molecule's FragmentationModel subgroups.
"""

import numpy as np

from rdkit import Chem

# from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


# def check_has_molecular_weight_right(
#     mol_object: Chem.rdchem.Mol,
#     mol_subgroups: dict,
#     model: FragmentationModel,
# ) -> bool:
#     """Check the molecular weight of the molecule using its functional groups.

#     Compares the RDKit molecular weight of the molecule to the computed
#     molecular weight from the functional groups. Returns True if both molecular
#     weights are equal with 0.5 u (half hydrogen atom) as atol of
#     numpy.allclose(). Also, the method will check if the molecule has negative
#     occurrences on its functional groups, also returning False in that case.

#     Parameters
#     ----------
#     mol_object : Chem.rdchem.Mol
#         RDKit Chem object
#     mol_subgroups : dict
#         FragmentationModel subgroups of the mol_object
#     model: FragmentationModel
#         FragmentationModel object.

#     Returns
#     -------
#     bool
#         True if RDKit and ugropy molecular weight are equal with a tolerance.
#     """
#     # check for negative occurrences
#     if not all(occurrence > 0 for occurrence in mol_subgroups.values()):
#         return False

#     # rdkit molecular weight
#     rdkit_mw = Descriptors.MolWt(mol_object)

#     # Molecular weight from functional groups
#     mws = model.subgroups.loc[
#         list(mol_subgroups.keys()), "molecular_weight"
#     ].to_numpy()

#     func_group_mw = np.dot(mws, list(mol_subgroups.values()))

#     return np.allclose(rdkit_mw, func_group_mw, atol=0.5)


# def check_can_fit_atoms(
#     mol_object: Chem.rdchem.Mol,
#     mol_subgroups: dict,
#     model: FragmentationModel,
# ) -> bool:
#     """Check if a solution can be fitted in the mol_object atoms.

#     Parameters
#     ----------
#     mol_object : Chem.rdchem.Mol
#         RDKit Mol object.
#     mol_subgroups : dict
#         Subgroups of mol_object.
#     model: FragmentationModel
#         FragmentationModel object.

#     Returns
#     -------
#     bool
#         True if the solution can be fitted.
#     """
#     if fit_atoms(mol_object, mol_subgroups, model):
#         return True
#     else:
#         return False


def check_has_overlapping_groups(
    mol_object: Chem.rdchem.Mol,
    mol_subgroups: dict,
) -> tuple[bool, np.ndarray]:
    """Check if the groups detection overlapping groups.

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Mol object.
    mol_subgroups : dict
        Subgroups of mol_object with the atoms indexes of each detection.
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    tuple[bool, np.ndarray]
        True if the groups detection has overlapping groups and the indexes of
        the overlapped atoms.
    """
    n_atoms = mol_object.GetNumAtoms()

    # Count the number of times an atom is in a group
    atoms = np.zeros(n_atoms)

    for indexes in mol_subgroups.values():
        np.add.at(atoms, np.array(indexes).flatten(), 1)

    overlapped_atoms = np.argwhere(atoms > 1).flatten()

    # Separate the fragments that participate in the overlapped atoms
    fragments = {}

    for oatom in overlapped_atoms:
        ...


    if np.size(overlapped_atoms) > 0:
        return True, overlapped_atoms
    else:
        return False, np.array([])
