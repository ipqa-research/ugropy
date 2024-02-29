"""detect_groups module."""

from typing import List

from rdkit import Chem

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


def detect_groups(
    mol_object: Chem.rdchem.Mol, model: FragmentationModel
) -> tuple[List[str], List[int]]:
    """Detect present functional groups in the mol_object molecule.

    Asks for each functional group in the subgroups DataFrame of the
    FragmentationModel using the SMARTS representation of the functional group.
    Then, returns the detected groups and the number of occurrences.

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Mol object
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    tuple[List[str], List[int]]
        Functional groups, functional group occurrences
    """
    groups = []
    occurrences = []

    for group in model.subgroups.index:
        matches = group_matches(mol_object, group, model)
        how_many_matches = len(matches)

        if how_many_matches > 0:
            groups.append(group)
            occurrences.append(int(how_many_matches))

    return groups, occurrences


def group_matches(
    mol_object: Chem.rdchem.Mol,
    group: str,
    model: FragmentationModel,
    action: str = "detection",
) -> tuple:
    """Obtain the group matches in mol_object.

    Given a functional group (group), a FragmentationModel (subgroups) and a
    RDKit Mol object (mol_object), obtains the SubstructMatches in chem_object
    returns a tuple of tuples containing the atoms that participate in the
    "group" substructure (return of the RDKit GetSubstructMatches function).
    The length of the return tuple is equal to the number of matches of the
    group in the mol_object.

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Mol object
    group : str
        String of the subgroup. E.g: 'CH3'
    model: FragmentationModel
        FragmentationModel object.
    action: str, optional
        Two options are possible: 'detection' or 'fit'. Choose the SMARTS
        representation according to the task.

    Returns
    -------
    tuple
        Return of the RDKit GetSubstructMatches function.
    """
    if action == "detection":
        mols = model.detection_mols[group]
    elif action == "fit":
        mols = model.fit_mols[group]
    else:
        raise ValueError(f"{action} not valid, use 'detection' or 'fit'")

    for mol in mols:
        matches = mol_object.GetSubstructMatches(mol)

        if len(matches) != 0:
            return matches

    return matches
