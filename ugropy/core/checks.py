"""check module.

The module contains the necessary checks to corroborate the success of the
algorithm to obtain the molecule's UNIFAC subgroups.

check_has_molecular_weight_right: check the molecular weight of the molecule.

check_has_composed: check if the molecule has composed structures.

check_has_hidden_ch2_ch: check if the molecule has CH2 or CH hidden in composed
structures.
"""

import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel

from .detect_model_groups import group_matches


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
) -> bool:
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
    composed_structures = model.subgroups[
        model.subgroups["composed"] == "y"
    ].index

    return any(composed in mol_subgroups for composed in composed_structures)


def check_has_hiden(
    mol_object: Chem.rdchem.Mol,
    mol_subgroups: dict,
    model: FragmentationModel,
) -> bool:
    """Check for hidden subgroups in composed structures.

    The principal subgroups that can be hidden in composed structures for the
    models UNIFAC, PSRK and Dortmund are CH2 and CH2. The algorithm checks that
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
        # import ipdb; ipdb.set_trace(cond=(candidate=="CH2"))
        misscount = 0

        exposed_candidates = mol_subgroups.get(candidate, 0)

        all_candidates_atoms = group_matches(mol_object, candidate, model)
        all_candidates_atoms = np.array(all_candidates_atoms).flatten()

        hideouts_atoms = np.array([])
        for hideout in model.hideouts.loc[candidate].values.flatten():
            if hideout in mol_subgroups.keys():
                atoms = group_matches(mol_object, hideout, model)

                # TODO: make documentation about the next if
                if len(atoms) > mol_subgroups[hideout] and model.subgroups.loc[hideout, "composed"] == "n":
                    misscount += len(atoms) - mol_subgroups[hideout]

                atoms = np.array(atoms).flatten()
                hideouts_atoms = np.append(hideouts_atoms, atoms)

        candidate_diff = np.setdiff1d(all_candidates_atoms, hideouts_atoms)

        if len(candidate_diff) + misscount != exposed_candidates:
            return True

    return False
