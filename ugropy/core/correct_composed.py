import json

from itertools import combinations

from rdkit import Chem

import numpy as np

import pandas as pd

from .checks import check_molecular_weight


def correct_composed(
    chem_object: Chem.rdchem.Mol, 
    chem_subgroups: dict, 
    subgroups: pd.DataFrame
):
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
            groups_to_decompose=combination
        )

        # Did the correction work?
        check_mw = check_molecular_weight(
            chem_object=chem_object,
            chem_subgroups=correction,
            subgroups=subgroups
        )
        
        if check_mw:
            successfull_corrections = np.append(
                successfull_corrections, 
                correction
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
        # TODO: Check this...
        subgroups_used = np.array(
            [np.sum(list(d.values())) for d in unique_corrections]
        )
        min_subgroups_used = np.min(subgroups_used)
        idx_min_lens = np.where(subgroups_used == min_subgroups_used)[0]
        dicts_with_min_len = unique_corrections[idx_min_lens]

        return dicts_with_min_len


def apply_decompose_correction(
        chem_subgroups: dict, 
        subgroups: pd.DataFrame,
        groups_to_decompose: tuple
    ):
    """

    Parameters
    ----------
    chem_subgroups : dict
        _description_
    subgroups : pd.DataFrame
        _description_
    groups_to_decompose : tuple[str]
        _description_

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
