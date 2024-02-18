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
    ACCH2 can be decomposed in the AC and CH2 groups.

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


def check_has_hidden_ch2_ch(
    mol_object: Chem.rdchem.Mol,
    mol_subgroups: dict,
    model: FragmentationModel,
) -> bool:
    """Check for hidden CH2 and CH subgroups in composed structures.

    The algorithm checks that the number of CH2 and CH groups in mol_subgroups
    dictionary is equal to the number of free CH2 and CH. If these numbers are
    not equal reveals that the CH2 and CH are hidden in composed structures,
    eg: ACCH2, ACCH. This phenomenon occurs when two subgroups fight for the
    same CH2 or CH. For example the molecule:

    CCCC1=CC=C(COC(C)(C)C)C=C1

    Here an ACCH2 and a CH2O are fighting to have the same CH2. But since there
    is a free CH2 in the molecule, the molecule prefers to keep both ACCH2 and
    CH2O groups without any free CH2 subgroup. This check searches for all CH2
    and counts all the CH2 that are participating in a CH2 hideout (ACCH2 and
    CH2O are examples of hideouts). The algorithm notices that there is one
    free CH2 and there are zero free CH2 groups in the mol_subgroups dictionary
    and returns 'True' (it has a hidden CH2 or CH).

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
        True if has hidden CH2 or CH subgroups.
    """
    ch2_num = mol_subgroups.get("CH2", 0)
    ch_num = mol_subgroups.get("CH", 0)

    all_ch2_atoms = group_matches(mol_object, "CH2", model)
    all_ch2_atoms = np.array(all_ch2_atoms).flatten()

    all_ch_atoms = group_matches(mol_object, "CH", model)
    all_ch_atoms = np.array(all_ch_atoms).flatten()

    ch2_hidings_atoms = np.array([])
    for hideout in model.ch2_hideouts:
        if hideout in mol_subgroups.keys():
            hidings = group_matches(mol_object, hideout, model)
            hidings = np.array(hidings).flatten()
            ch2_hidings_atoms = np.append(ch2_hidings_atoms, hidings)

    ch_hidings_atoms = np.array([])
    for hideout in model.ch_hideouts:
        if hideout in mol_subgroups.keys():
            hidings = group_matches(mol_object, hideout, model)
            hidings = np.array(hidings).flatten()
            ch_hidings_atoms = np.append(ch_hidings_atoms, hidings)

    ch2_diff = np.setdiff1d(all_ch2_atoms, ch2_hidings_atoms)
    ch_diff = np.setdiff1d(all_ch_atoms, ch_hidings_atoms)

    if len(ch2_diff) != ch2_num:
        return True
    elif len(ch_diff) != ch_num:
        return True
    else:
        return False
