import json

from itertools import combinations

from rdkit import Chem
from rdkit.Chem import Descriptors

import numpy as np
from numpy.typing import NDArray

import pandas as pd


def get_groups(
    chem_object: Chem.rdchem.Mol, 
    subgroups: pd.DataFrame, 
    subgroups_matrix: pd.DataFrame, 
    problematic_structures: pd.DataFrame
):
    # Shorter names for DataFrames:
    df = subgroups.copy()
    dfm = subgroups_matrix.copy()
    dfp = problematic_structures.copy()

    # Groups found and the occurrence number of each one in chem_object.
    groups, groups_ocurrences = detect_groups(
        chem_object=chem_object, 
        subgroups=df
    )

    # Filters the subgroups matrix into the detected groups only.
    dff = dfm.loc[groups][groups]
    
    # Multiply each group row by the occurrences of that group.
    dff = dff.mul(groups_ocurrences, axis=0)

    # Correction of problematic structures in dff.
    dff_corrected = correct_problematics(
        chem_object=chem_object,
        filtered_subgroups=dff,
        problematic_structures=dfp,
    )

    # Calculate the number of each functional group.
    dff_sum = dff_corrected.sum(axis=0)
    dff_sum.replace(0, pd.NA, inplace=True)
    chem_subgroups = dff_sum.dropna()
    chem_subgroups = chem_subgroups.to_dict()

    if chem_subgroups == {}:
        # No functional groups detected for the molecule. Example: hydrogen
        # peroxide.
        return chem_subgroups
    
    # Check for composed structures.
    if check_molecular_weight(chem_object=chem_object, chem_subgroups=chem_subgroups, subgroups=df):
        return chem_subgroups
    else:
        chem_subgroups = correct_composed(
            chem_object=chem_object,
            molecule_func_groups=chem_subgroups,
            subgroups=df
        )
        return chem_subgroups


def detect_groups(
    chem_object: Chem.rdchem.Mol, 
    subgroups: pd.DataFrame
) -> tuple[NDArray, NDArray]:
    """Detect present functional groups in the chem_object molecule.
    
    Asks for each functional group in the subgroups DataFrame using the SMARTS 
    representation of the functional group. Then, stores the detected group and
    the number of occurrences.

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Chem object
    subgroups : pd.DataFrame
        DataFrame with the UNIFAC model's subgroups

    Returns
    -------
    tuple[NDArray, NDArray]
        
    """
    df = subgroups.copy()
    
    groups = np.array([])
    groups_ocurrences = np.array([])
    
    for group in df.index:
        smarts = df.loc[group]["smarts"]

        func_group = Chem.MolFromSmarts(smarts)
        matches = chem_object.GetSubstructMatches(func_group)
        how_many_matches = len(matches)
        
        if how_many_matches > 0:
            groups = np.append(groups, group)
            groups_ocurrences = np.append(
                groups_ocurrences, 
                how_many_matches
            ).astype(int)
            
    return groups, groups_ocurrences


def correct_problematics(
    chem_object: Chem.rdchem.Mol, 
    filtered_subgroups:pd.DataFrame, 
    problematic_structures: pd.DataFrame
) -> pd.DataFrame:
    """Corrects problematic structures in chem_object.
    
    The algorithm of the function detect_groups may have some troubles with
    some chemical structures. For example, in an ester group, it will detect
    also an ether group. This problem can be handled by changing the SMARTS
    representation of the ether group to something like:
    
    [CH3][O][[#6]&!$([C]=O)]
    
    With this SMARTS representation of the CH3O UNIFAC group, we are specifying
    that this functional group it's a methyl group (CH3) bounded to oxygen by a
    simple covalent bond, and that oxygen is bonded to any carbon but a carbon 
    that it's double bonded to an oxygen. This should avoid detecting an ether 
    group in an ester group. But, consider the structure of the molecule ethyl 
    methyl carbonate (PubChem CID 522046). That molecule has both an ester 
    group and an ether group, and the previous smarts representation will not 
    detect the ether group that we want to be detected. This problem defines 
    what a problematic structure is for the ugropy library. Maybe it is a 
    SMARTS representation that well behaves in these situations, but it's 
    easier to make a list of these problematic structures and correct them with
    a function than generate a complex SMARTS representation of a functional 
    group.

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Chem object
    filtered_subgroups : pd.DataFrame
        Subgroups matrix DataFrame filtered for groups present in chem_object.
    problematic_structures : pd.DataFrame
        Problematic structures DataFrame.

    Returns
    -------
    pd.DataFrame
        Subgroups matrix DataFrame filtered with corrected problematic 
        structures. 
    """
    dff = filtered_subgroups.copy()
    dfp = problematic_structures.copy()
    
    for smarts in dfp.index:
        structure = Chem.MolFromSmarts(smarts)
        matches = chem_object.GetSubstructMatches(structure)
        how_many_problems = len(matches)

        if how_many_problems > 0:
            problm_dict = json.loads(dfp.loc[smarts].contribute)

            for grp in problm_dict.keys():
                dff.loc[grp][grp] += problm_dict[grp] * how_many_problems
                
    return dff


def check_molecular_weight(
    chem_object: Chem.rdchem.Mol, 
    chem_subgroups: dict,
    subgroups: pd.DataFrame,
) -> bool:
    """Check the molecular weight of the molecule using its functional groups.
    
    Compares the RDKit molecular weight of the molecule to the computed
    colecular weight from the functional groups. Returns True if both molecular
    weight are equals with 0.5 u (half hydrogen atom) as atol of 
    numpy.allclose().

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Chem object
    chem_subgroups : dict
        dict with the UNIFAC subgroups of the chem_object
    subgroups : pd.DataFrame
        DataFrame with the UNIFAC model's subgroups
    tolerance : float
        Allowed difference between RDKit and ugropy molecular weights

    Returns
    -------
    bool
        True if RDKit and ugropy molecular weight are equals with a tolerance.
    """
    dff = chem_subgroups.copy()

    # check for negative occurrences
    if not all(occurrence > 0 for occurrence in dff.values()):
        return False
    
    # rdkit molecular weight
    rdkit_mw = Descriptors.MolWt(chem_object)
    
    # Get molecular weight from functional groups
    func_group_mw = 0

    tolerance = 0.5 # 1/2 hydrogen atom.
    
    for group in chem_subgroups.keys():
        func_group_mw += subgroups.loc[group]["molecular_weight"] * chem_subgroups[group]
    
    return np.allclose(rdkit_mw, func_group_mw, atol=tolerance)


def correct_composed(chem_object: Chem.rdchem.Mol, molecule_func_groups: dict, subgroups: pd.DataFrame):
    chem_subgroups = molecule_func_groups.copy()
    df = subgroups.copy()

    # Get composed structures from subgroups DataFrame
    composed_structures = df[df["composed"] == "y"].index

    # Create list of composed structures present in molecule's func_g. If a 
    # composed structure has occurrence 2, will appear twice in the list
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

    if len(combinatory_list) == 0:
        # There is no possible correction and probably no way to represent the
        # molecule by UNIFAC's functional groups.
        return {}


    # Try by brute force all combinatories and store the successfull ones
    successfull_corrections = np.array([])

    for combination in combinatory_list:
        # Get subgroups with decomposed structures
        correction = apply_decompose_correction(
            chem_subgroups=chem_subgroups,
            subgroups=df,
            combination=combination
        )

        # Did the correction work?
        did_it_work = check_molecular_weight(
            chem_object=chem_object,
            chem_subgroups=correction,
            subgroups=subgroups
        )

        if did_it_work:
            successfull_corrections = np.append(successfull_corrections, correction)

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

        return dicts_with_min_len


def check_ghosts(
    chem_object: Chem.rdchem.Mol, 
    molecule_func_groups: dict, 
    subgroups: pd.DataFrame
) -> bool:
    func_groups = molecule_func_groups.copy()
    df = subgroups.copy()
    
    # Get composed structures from subgroups DataFrame
    composed_structures = df[df["composed"] == "y"].index
    
    # Get composed structures in molecule
    composed_in_mol = np.array([])
    for composed in composed_structures:
        if composed in func_groups.key():
            composed_in_mol = np.append(composed_in_mol, composed)
            
    # Get uncomposed structures in molecule
    not_composed_in_mol = np.array([])
    for structure in func_groups:
        if structure not in composed_structures:
            not_composed_in_mol = np.append(not_composed_in_mol, structure)
            
    # Get atoms participating in composed structures
    composed_atoms = {}
    for fgrp in composed_in_mol:
        smarts = df.loc[fgrp]["smarts"]
        mol = Chem.MolFromSmarts(smarts)
        matches = chem_object.GetSubstructMatches(mol)
        composed_atoms.update({fgrp: matches})
        
    # Get atoms participatins in no composed structures
    no_composed_atoms = {}
    for fgrp in not_composed_in_mol:
        smarts = df.loc[fgrp]["smarts"]
        mol = Chem.MolFromSmarts(fgrp)
        matches = chem_object.GetSubstructMatches(mol)
        no_composed_atoms.update({fgrp: matches})
        
    # Get supect atoms
    for c_fgrp in composed_in_mol:
        # Get decomposed funcitonal groups
        contribute = json.loads(df.loc[c_fgrp].contribute)
        decomposed = contribute.pop(c_fgrp)
        
        decomposed_atoms = {}
        
        # Need the decomposed atoms in structure
        for dcmp in decomposed.keys():
            smarts = df.loc[dcmp]["smarts"]
            mol = Chem.MolFromSmarts(smarts)
            matches = chem_object.GetSubstructMatches(mol)
            decomposed_atoms.update({dcmp: matches})


def apply_decompose_correction(
        chem_subgroups: dict, 
        subgroups: pd.DataFrame,
        combination: tuple
    ): 
    
    chm_grps = chem_subgroups.copy()
    df = subgroups.copy()

    for group in combination:
        contribute_dict = json.loads(df.loc[group].contribute)
        for grp in contribute_dict.keys():
            try:
                chm_grps[grp] += -1 * contribute_dict[grp]
            except KeyError:
                chm_grps[grp] = -1 * contribute_dict[grp]

    # Eliminate occurrences == 0
    groups_corrected = {key: value for key, value in chm_grps.items() if value > 0}

    return groups_corrected
