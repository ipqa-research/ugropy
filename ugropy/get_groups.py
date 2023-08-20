import json

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
    dff_final = dff_sum.dropna()
    
    # Check for composed structures.
    if check_molecular_weight(chem_object=chem_object, chem_object_subgroups=dff_final, subgroups=df):
        return dff_final.to_dict()
    #return dff_final.to_dict()


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
    chem_object_subgroups: pd.DataFrame,
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
    chem_object_subgroups : pd.DataFrame
        DataFrame with the UNIFAC subgroups of the chem_object
    subgroups : pd.DataFrame
        DataFrame with the UNIFAC model's subgroups
    tolerance : float
        Allowed difference between RDKit and ugropy molecular weights

    Returns
    -------
    bool
        True if RDKit and ugropy molecular weight are equals with a tolerance.
    """
    df = subgroups.copy()
    dff = chem_object_subgroups.copy()
    
    # rdkit molecular weight
    rdkit_mw = Descriptors.MolWt(chem_object)
    
    # Get molecular weight from functional groups
    func_group_mw = 0

    tolerance = 0.5 # 1/2 hydrogen atom.
    
    for group in dff.index:
        func_group_mw += df.loc[group]["molecular_weight"] * dff[group]
    
    return np.allclose(rdkit_mw, func_group_mw, atol=tolerance)
    