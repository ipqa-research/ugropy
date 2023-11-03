"""detect_groups module."""
import numpy as np
from numpy.typing import NDArray

import pandas as pd

from rdkit import Chem


def detect_groups(
    chem_object: Chem.rdchem.Mol, subgroups: pd.DataFrame
) -> tuple[NDArray, NDArray]:
    """Detect present functional groups in the chem_object molecule.

    Asks for each functional group in the subgroups DataFrame using the SMARTS
    representation of the functional group. Then, returns the detected groups
    and the number of occurrences.

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Mol object
    subgroups : pd.DataFrame
        DataFrame with the UNIFAC model's subgroups

    Returns
    -------
    tuple[NDArray, NDArray]
        Functional groups, functional group occurrences
    """
    groups = np.array([])
    occurrences = np.array([])

    for group in subgroups.index:
        matches = group_matches(chem_object, group, subgroups)
        how_many_matches = len(matches)

        if how_many_matches > 0:
            groups = np.append(groups, group)
            occurrences = np.append(occurrences, how_many_matches).astype(int)

    return groups, occurrences


def group_matches(
    chem_object: Chem.rdchem.Mol, group: str, subgroups: pd.DataFrame
) -> tuple:
    """Obtain the group matches in chem_object.

    Given a functional group (group), a subgroup DataFrame (subgroups) and a
    RDKit Mol object (chem_object), obtains the SubstructMatches in chem_object
    returns a tuple of tuples containing the atoms that participate in the
    "group" substructure (return of the RDKit GetSubstructMatches function).
    The length of the return tuple is equal to the number of matches of the
    group in the chem_object.

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Mol object
    group : str
        String of the subgroup. E.g: 'CH3'
    subgroups : pd.DataFrame
        DataFrame with the UNIFAC model's subgroups

    Returns
    -------
    tuple
        Return of the RDKit GetSubstructMatches function.
    """
    smarts = subgroups.loc[group, "smarts"]
    func_group = Chem.MolFromSmarts(smarts)

    matches = chem_object.GetSubstructMatches(func_group)

    return matches
