"""correct_composed module."""
import json
from itertools import combinations

import numpy as np

import pandas as pd

from rdkit import Chem

from .checks import check_has_hidden_ch2_ch, check_has_molecular_weight_right


def correct_composed(
    chem_object: Chem.rdchem.Mol,
    chem_subgroups: dict,
    subgroups: pd.DataFrame,
    ch2_hideouts: pd.DataFrame,
    ch_hideouts: pd.DataFrame,
) -> dict:
    """Correct composed structures.

    A priori is not easy to recognize what composed structures in
    chem_subgroups need to be decomposed to correct the solution. By that, all
    the combinations are tried. For example, a molecule that can't be solved
    has one ACCH2 and two ACCH composed structures. The decomposition
    combinatory will be:

    [[ACCH2], [ACCH], [ACCH2, ACCH], [ACCH, ACCH], [ACCH2, ACCH, ACCH]]

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Mol object.
    chem_subgroups : dict
        Molecule's UNIFAC subgroups.
    subgroups : pd.DataFrame
        DataFrame of all UNIFAC's subgroups.
    ch2_hideouts : pandas.DataFrame
        DataFrame of all posible CH2 group hidings.
    ch_hideouts : pandas.DataFrame
        DataFrame of all posible CH group hidings.

    Returns
    -------
    dict or list[dict]
        Corrected subgroups due to decomposing composed structures.
    """
    # Create list of composed structures present in molecule's func_g. If a
    # composed structure has occurrence 2, will appear twice in the list
    composed_structures = subgroups[subgroups["composed"] == "y"].index
    comp_in_mol = np.array([])

    for stru in composed_structures:
        if stru in chem_subgroups.keys():
            comp_in_mol = np.append(comp_in_mol, [stru] * chem_subgroups[stru])

    # Create a combinatory list of the possible decomposition of the structures
    combinatory_list = []
    for i in range(1, len(comp_in_mol) + 1):
        # turn into set eliminates duplicated corrections combinations.
        for combinatory in set(combinations(comp_in_mol, i)):
            combinatory_list.append(combinatory)

    # Try by brute force all combinatories and store the successfull ones
    successfull_corrections = np.array([])

    for combination in combinatory_list:
        # Get subgroups with decomposed structures
        correction = apply_decompose_correction(
            chem_subgroups=chem_subgroups,
            subgroups=subgroups,
            groups_to_decompose=combination,
        )

        # Did the correction work?
        right_mw = check_has_molecular_weight_right(
            chem_object=chem_object,
            chem_subgroups=correction,
            subgroups=subgroups,
        )

        has_hidden = check_has_hidden_ch2_ch(
            chem_object=chem_object,
            chem_subgroups=correction,
            subgroups=subgroups,
            ch2_hideouts=ch2_hideouts,
            ch_hideouts=ch_hideouts,
        )

        if right_mw and not has_hidden:
            successfull_corrections = np.append(
                successfull_corrections, correction
            )

    # No posible correction found, can't represent molecule with func groups
    if len(successfull_corrections) == 0:
        return {}

    # Get rid of duplicated successfull_corrections
    dict_strs = np.array([str(d) for d in successfull_corrections])
    unique_indices = np.unique(dict_strs, return_index=True)[1]
    unique_corrections = successfull_corrections[unique_indices]

    # Return decomposed subgroups solutions
    if len(unique_corrections) == 1:
        # Unique solution found
        return unique_corrections[0]
    else:
        # Find the solution/s that uses the minimun number of functional groups
        subgroups_used = np.array(
            [np.sum(list(d.values())) for d in unique_corrections]
        )
        min_subgroups_used = np.min(subgroups_used)
        idx_min_lens = np.where(subgroups_used == min_subgroups_used)[0]
        dicts_with_min_len = unique_corrections[idx_min_lens]

        if len(dicts_with_min_len) == 1:
            # Solutions number turn into one
            return dicts_with_min_len[0]
        else:
            return dicts_with_min_len


def apply_decompose_correction(
    chem_subgroups: dict, subgroups: pd.DataFrame, groups_to_decompose: tuple
) -> dict:
    """Decompose composed structures in chem_subgroups.

    The function receives a tuple of groups to decompose and applies the
    corresponding correction. For example, if the function receives the order
    of decomposing an ACCH2, the function subtracts an ACCH2 group from
    chem_subgroups, and then adds an AC and a CH2 group.

    Parameters
    ----------
    chem_subgroups : dict
        Molecule's UNIFAC subgroups.
    subgroups : pd.DataFrame
        DataFrame with the UNIFAC's model subgroups.
    groups_to_decompose : tuple[str]
        Tuple with all the composed structures to decompose.

    Returns
    -------
    dict
        Functional groups dictionary with decomposed structures.
    """
    chm_grps = chem_subgroups.copy()

    for group in groups_to_decompose:
        contribute_dict = json.loads(subgroups.loc[group].contribute)
        for grp in contribute_dict.keys():
            try:
                chm_grps[grp] += -1 * contribute_dict[grp]
            except KeyError:
                chm_grps[grp] = -1 * contribute_dict[grp]

    # Eliminate occurrences == 0
    groups_corrected = {
        key: value for key, value in chm_grps.items() if value != 0
    }

    return groups_corrected
