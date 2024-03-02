"""correct_composed module."""

import json
from itertools import combinations

import numpy as np

from rdkit import Chem

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel

from .checks import (
    check_can_fit_atoms,
    check_has_composed_overlapping,
    check_has_hiden,
    check_has_molecular_weight_right,
)


def correct_composed(
    mol_object: Chem.rdchem.Mol,
    mol_subgroups: dict,
    model: FragmentationModel,
) -> dict:
    """Correct composed structures.

    A priori is not easy to recognize what composed structures in
    mol_subgroups need to be decomposed to correct the solution. By that, all
    the combinations are tried. For example, a molecule that can't be solved
    has one ACCH2 and two ACCH composed structures. The decomposition
    combinatory will be:

    [[ACCH2], [ACCH], [ACCH2, ACCH], [ACCH, ACCH], [ACCH2, ACCH, ACCH]]

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Mol object.
    mol_subgroups : dict
        Molecule's FragmentationModel subgroups.
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    dict or list[dict]
        Corrected subgroups due to decomposing composed structures.
    """
    # =========================================================================
    # A list with all the composed structures present in mol_subgroups
    # =========================================================================
    composed_structures = [
        stru
        for stru in mol_subgroups.keys()
        if model.subgroups.loc[stru, "composed"] == "y"
    ]

    # =========================================================================
    # Creates a list with the composed structures in mol_subgroups but each
    # composed structure is repetead a number of times equals to the occurences
    # In te molecule. For example in UNIFAC {"ACCH2": 3, "CH3": 1, "ACCH": 1}
    # should generate:
    #
    # ["ACCH2", "ACCH2", "ACCH2", "ACCH"]
    # =========================================================================
    composed_in_mol = np.array(
        [
            stru
            for stru in composed_structures
            for _ in range(mol_subgroups[stru])
        ]
    )

    # =========================================================================
    # Create the combinatory list as explainde in the funcion documentation.
    # =========================================================================
    combinatory_list = []
    for i in range(1, len(composed_in_mol) + 1):
        combinatory_list.extend(set(combinations(composed_in_mol, i)))

    # =========================================================================
    # Try by brute force all combinatories and store the successfull ones
    # =========================================================================
    successfull_corrections = []

    for combination in combinatory_list:
        # Get subgroups with decomposed structures
        correction = apply_decompose_correction(
            mol_subgroups=mol_subgroups,
            groups_to_decompose=combination,
            model=model,
        )

        # Did the correction work?
        right_mw = check_has_molecular_weight_right(
            mol_object=mol_object,
            mol_subgroups=correction,
            model=model,
        )

        if not right_mw:
            continue

        has_overlap = check_has_composed_overlapping(
            mol_object, correction, model
        )

        if has_overlap:
            continue

        has_hiden = check_has_hiden(mol_object, correction, model)

        if has_hiden:
            continue

        can_fit = check_can_fit_atoms(
            mol_object=mol_object,
            mol_subgroups=correction,
            model=model,
        )

        if can_fit:
            successfull_corrections.append(correction)

    successfull_corrections = np.array(successfull_corrections)

    # No posible correction found, can't represent molecule with func groups
    if len(successfull_corrections) == 0:
        return {}

    # =========================================================================
    # Get rid of duplicated successfull_corrections
    # =========================================================================
    dict_strs = np.array([str(d) for d in successfull_corrections])
    unique_indices = np.unique(dict_strs, return_index=True)[1]
    unique_corrections = successfull_corrections[unique_indices]

    # =========================================================================
    # Return decomposed subgroup solution
    # =========================================================================
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
    mol_subgroups: dict,
    groups_to_decompose: tuple,
    model: FragmentationModel,
) -> dict:
    """Decompose composed structures in mol_subgroups.

    The function receives a tuple of groups to decompose and applies the
    corresponding correction. For example, if the function receives the order
    of decomposing an ACCH2, the function subtracts an ACCH2 group from
    mol_subgroups, and then adds an AC and a CH2 group.

    Parameters
    ----------
    mol_subgroups : dict
        Molecule's FragmentationModel subgroups.
    groups_to_decompose : tuple[str]
        Tuple with all the composed structures to decompose.
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    dict
        Functional groups dictionary with decomposed structures.
    """
    chm_grps = mol_subgroups.copy()

    for group in groups_to_decompose:
        contribute_dict = json.loads(model.subgroups.loc[group, "contribute"])
        for grp, contribution in contribute_dict.items():
            chm_grps[grp] = chm_grps.get(grp, 0) - 1 * contribution

    # Eliminate occurrences == 0
    groups_corrected = {
        key: value for key, value in chm_grps.items() if value != 0
    }

    return groups_corrected
