"""check module.

The module contains the necessary checks to corroborate the success of the
algorithm to obtain the molecule's UNIFAC subgroups.

check_has_molecular_weight_right: check the molecular weight of the molecule.

check_has_composed: check if the molecule has composed structures.

check_has_hidden_ch2_ch: check if the molecule has CH2 or CH hidden in composed
structures.
"""
import numpy as np

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors

from .detect_groups import group_matches


def check_has_molecular_weight_right(
    chem_object: Chem.rdchem.Mol,
    chem_subgroups: dict,
    subgroups: pd.DataFrame,
) -> bool:
    """Check the molecular weight of the molecule using its functional groups.

    Compares the RDKit molecular weight of the molecule to the computed
    molecular weight from the functional groups. Returns True if both molecular
    weights are equal with 0.5 u (half hydrogen atom) as atol of
    numpy.allclose(). Also, the method will check if the molecule has negative
    occurrences on its functional groups, also returning False in that case.

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Chem object
    chem_subgroups : dict
        UNIFAC subgroups of the chem_object
    subgroups : pd.DataFrame
        DataFrame with the UNIFAC model's subgroups
    tolerance : float
        Allowed difference between RDKit and ugropy molecular weights

    Returns
    -------
    bool
        True if RDKit and ugropy molecular weight are equal with a tolerance.
    """
    # check for negative occurrences
    if not all(occurrence > 0 for occurrence in chem_subgroups.values()):
        return False

    # rdkit molecular weight
    rdkit_mw = Descriptors.MolWt(chem_object)

    # Molecular weight from functional groups
    mws = subgroups.loc[chem_subgroups.keys()]["molecular_weight"].to_numpy()
    func_group_mw = np.dot(mws, list(chem_subgroups.values()))

    return np.allclose(rdkit_mw, func_group_mw, atol=0.5)


def check_has_composed(
    chem_subgroups: dict,
    subgroups: pd.DataFrame,
) -> bool:
    """Check if the molecule has composed structures.

    A composed structure is a subgroup of UNIFAC that can be decomposed into
    two or more UNIFAC subgroups. For example, ACCH2 can be decomposed in the
    AC and CH2 groups.

    Parameters
    ----------
    chem_subgroups : dict
        Dictionary with the detected subgroups.
    subgroups : pd.DataFrame
        Complete DataFrame of UNIFAC's subgroups.

    Returns
    -------
    bool
        True if the molecule has composed structures.
    """
    composed_structures = subgroups[subgroups["composed"] == "y"].index

    for composed in composed_structures:
        if composed in chem_subgroups.keys():
            return True

    return False


def check_has_hidden_ch2_ch(
    chem_object: Chem.rdchem.Mol,
    chem_subgroups: dict,
    subgroups: pd.DataFrame,
    ch2_hideouts: pd.DataFrame,
    ch_hideouts: pd.DataFrame,
) -> bool:
    """Check for hidden CH2 and CH subgroups in composed structures.

    The algorithm checks that the number of CH2 and CH groups in chem_subgroups
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
    free CH2 and there are zero free CH2 groups in the chem_subgroups
    dictionary and returns 'True' (it has a hidden CH2 or CH).

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Mol object.
    chem_subgroups : dict
        Subgroups of chem_object.
    subgroups : pd.DataFrame
        Complete DataFrame of UNIFAC's subgroups.
    ch2_hideouts : pandas.DataFrame
        DataFrame of all posible CH2 group hidings.
    ch_hideouts : pandas.DataFrame
        DataFrame of all posible CH group hidings.

    Returns
    -------
    bool
        True if has hidden CH2 or CH subgroups.
    """
    try:
        ch2_num = chem_subgroups["CH2"]
    except KeyError:
        ch2_num = 0

    try:
        ch_num = chem_subgroups["CH"]
    except KeyError:
        ch_num = 0

    all_ch2_atoms = group_matches(chem_object, "CH2", subgroups)
    all_ch2_atoms = np.array(all_ch2_atoms).flatten()

    all_ch_atoms = group_matches(chem_object, "CH", subgroups)
    all_ch_atoms = np.array(all_ch_atoms).flatten()

    ch2_hidings_atoms = np.array([])
    for hideout in ch2_hideouts:
        if hideout in chem_subgroups.keys():
            hidings = group_matches(chem_object, hideout, subgroups)
            hidings = np.array(hidings).flatten()
            ch2_hidings_atoms = np.append(ch2_hidings_atoms, hidings)

    ch_hidings_atoms = np.array([])
    for hideout in ch_hideouts:
        if hideout in chem_subgroups.keys():
            hidings = group_matches(chem_object, hideout, subgroups)
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
