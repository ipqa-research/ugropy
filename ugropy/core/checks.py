"""check module.

The module contains the necessary checks to corroborate the success of the
algorithm to obtain the molecule's FragmentationModel subgroups.
"""

import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel

from .detect_model_groups import group_matches
from .fit_atoms_indexes import fit_atoms


def check_has_molecular_weight_right(
    mol_object: Chem.rdchem.Mol,
    mol_subgroups: dict,
    model: FragmentationModel,
) -> bool:
    """Check the molecular weight of the molecule using its functional groups.

    Compares the RDKit molecular weight of the molecule to the computed
    molecular weight from the functional groups. Returns True if both molecular
    weights are equal with 0.5 u (half hydrogen atom) as atol of
    numpy.allclose(). Also, the method will check if the molecule has negative
    occurrences on its functional groups, also returning False in that case.

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Chem object
    mol_subgroups : dict
        FragmentationModel subgroups of the mol_object
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    bool
        True if RDKit and ugropy molecular weight are equal with a tolerance.
    """
    # check for negative occurrences
    if not all(occurrence > 0 for occurrence in mol_subgroups.values()):
        return False

    # rdkit molecular weight
    rdkit_mw = Descriptors.MolWt(mol_object)

    # Molecular weight from functional groups
    mws = model.subgroups.loc[
        list(mol_subgroups.keys()), "molecular_weight"
    ].to_numpy()

    func_group_mw = np.dot(mws, list(mol_subgroups.values()))

    return np.allclose(rdkit_mw, func_group_mw, atol=0.5)


def check_has_composed(
    mol_subgroups: dict,
    model: FragmentationModel,
) -> tuple[bool, np.ndarray]:
    """Check if the molecule has composed structures.

    A composed structure is a subgroup of FragmentationModel that can be
    decomposed into two or more FragmentationModel subgroups. For example,
    ACCH2 can be decomposed into the AC and CH2 groups.

    Parameters
    ----------
    mol_subgroups : dict
        Dictionary with the detected subgroups.
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    bool
        True if the molecule has composed structures.
    """
    composed_stru = model.subgroups[model.subgroups["composed"] == "y"].index
    composed_in_mol = np.intersect1d(composed_stru, list(mol_subgroups.keys()))
    return len(composed_in_mol) > 0, composed_in_mol


def check_has_hiden(
    mol_object: Chem.rdchem.Mol,
    mol_subgroups: dict,
    model: FragmentationModel,
) -> bool:
    """Check for hidden subgroups in composed structures.

    The principal subgroups that can be hidden in composed structures for the
    models UNIFAC, PSRK and Dortmund are CH2 and CH. The algorithm checks that
    the number of CH2 and CH groups in mol_subgroups dictionary is equal to the
    number of free CH2 and CH. If these numbers are not equal reveals that the
    CH2 and CH are hidden in composed structures, eg: ACCH2, ACCH. This
    phenomenon occurs when two subgroups fight for the same CH2 or CH. For
    example the molecule:

    CCCC1=CC=C(COC(C)(C)C)C=C1

    Here an ACCH2 and a CH2O are fighting to have the same CH2. But since there
    is a free CH2 in the molecule, the algorithm prefers to keep both ACCH2 and
    CH2O groups without any free CH2 subgroup. This check counts all the CH2
    that are participating in a CH2 hideout (ACCH2 and CH2O are examples of
    hideouts). The algorithm notices that there is one free CH2 and there are
    zero free CH2 groups in the mol_subgroups dictionary and returns 'True'
    (mol_object has a hidden CH2).

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Mol object.
    mol_subgroups : dict
        Subgroups of mol_object.
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    bool
        True if has hidden subgroups.
    """
    hiden_candidates = np.unique(model.hideouts.index.to_numpy())

    for candidate in hiden_candidates:
        exposed_candidates = mol_subgroups.get(candidate, 0)

        all_candidates_atoms = group_matches(mol_object, candidate, model)
        all_candidates_atoms = np.array(all_candidates_atoms).flatten()

        hideouts_atoms = np.array([])
        for hideout in model.hideouts.loc[candidate].values.flatten():
            if hideout in mol_subgroups.keys():
                atoms = group_matches(mol_object, hideout, model, "fit")

                atoms = np.array(atoms).flatten()
                hideouts_atoms = np.append(hideouts_atoms, atoms)

                much_matches = len(atoms) > mol_subgroups[hideout]
                no_comp = model.subgroups.loc[hideout, "composed"] == "n"

                # TODO: document this
                if much_matches and no_comp:
                    return False

        candidate_diff = np.setdiff1d(all_candidates_atoms, hideouts_atoms)

        if len(candidate_diff) != exposed_candidates:
            return True

    return False


def check_can_fit_atoms(
    mol_object: Chem.rdchem.Mol,
    mol_subgroups: dict,
    model: FragmentationModel,
) -> bool:
    """Check if a solution can be fitted in the mol_object atoms.

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Mol object.
    mol_subgroups : dict
        Subgroups of mol_object.
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    bool
        True if the solution can be fitted.
    """
    if fit_atoms(mol_object, mol_subgroups, model):
        return True
    else:
        return False


def check_has_composed_overlapping(
    mol_object: Chem.rdchem.Mol,
    mol_subgroups: dict,
    model: FragmentationModel,
) -> bool:
    """Check if in the solution are composed structures overlapping.

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Mol object.
    mol_subgroups : dict
        Subgroups of mol_object.
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    bool
        Treu if the solution has overlapping composed structures.
    """
    # =========================================================================
    # Count total number of composed in mol_subgroups
    # =========================================================================
    _, composed = check_has_composed(mol_subgroups=mol_subgroups, model=model)
    composed_in_subgroups = np.sum([mol_subgroups[gr] for gr in composed])

    # =========================================================================
    # Get atoms of composed and check overlaps
    # =========================================================================
    composed_atoms = [group_matches(mol_object, c, model) for c in composed]
    total_composed_matches = np.sum([len(c) for c in composed_atoms])

    overlapping_count = 0

    # Self overlapping
    for c_atoms in composed_atoms:
        atoms_array = np.array(c_atoms).flatten()
        _, counts = np.unique(atoms_array, return_counts=True)
        overlapping_count += np.sum(counts - 1)

    # Cross overlapping
    for i, i_atoms in enumerate(composed_atoms):
        for j_atoms in composed_atoms[i + 1 :]:  # noqa
            overlapping_count += np.sum(np.isin(i_atoms, j_atoms))

    response = composed_in_subgroups > (
        total_composed_matches - overlapping_count
    )

    return response
