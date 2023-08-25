from rdkit import Chem
from rdkit.Chem import Descriptors

import numpy as np

import pandas as pd

from .detect_groups import group_matches


def check_molecular_weight(
    chem_object: Chem.rdchem.Mol, 
    chem_subgroups: dict,
    subgroups: pd.DataFrame,
) -> bool:
    """Check the molecular weight of the molecule using its functional groups.
    
    Compares the RDKit molecular weight of the molecule to the computed 
    molecular weight from the functional groups. Returns True if both molecular
    weight are equals with 0.5 u (half hydrogen atom) as atol of 
    numpy.allclose(). Also the method will check if the molecule has negative
    occurrences on it's functional groups, also returning False in that case.

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
    # check for negative occurrences
    if not all(occurrence > 0 for occurrence in chem_subgroups.values()):
        return False
    
    # rdkit molecular weight
    rdkit_mw = Descriptors.MolWt(chem_object)
    
    # Molecular weight from functional groups
    mws = subgroups.loc[chem_subgroups.keys()]["molecular_weight"].to_numpy()
    func_group_mw = np.dot(mws, list(chem_subgroups.values()))
    
    return np.allclose(rdkit_mw, func_group_mw, atol=0.5)


def check_sneaky_ch2_ch(
    chem_object: Chem.rdchem.Mol, 
    chem_subgroups: dict, 
    subgroups: pd.DataFrame
) -> bool:
    subgroups_no_ch2_ch = chem_subgroups.copy()    
    
    ch2_matches = np.array(
        group_matches(chem_object, "CH2", subgroups)
    ).flatten()
    ch_matches = np.array(
        group_matches(chem_object, "CH", subgroups)
    ).flatten()

    # Get composed structures from subgroups DataFrame
    composed_structures = subgroups[subgroups["composed"] == "y"].index

    comp_in_mol = np.array([])
    for stru in composed_structures:
        if stru in chem_subgroups.keys():
            comp_in_mol = np.append(comp_in_mol, [stru] * chem_subgroups[stru])
            
    # Atoms of composed structures            
    grps_matches = []
    for group in chem_subgroups.keys():
        matches = group_matches(chem_object, group, subgroups)
        if len(matches) > 0:
            grps_matches.append(matches)
    
    # Flat composed matches
    f_composed_matches = np.array([atom for subtuple in grps_matches for atom in subtuple])
    
    lonely_ch2 = np.setdiff1d(ch2_matches, f_composed_matches)
    lonely_ch = np.setdiff1d(ch_matches, f_composed_matches)
    
    try:
        ch2_num = chem_subgroups["CH2"]
    except KeyError:
        ch2_num = 0
        
    try:
        ch_num = chem_subgroups["CH"]
    except KeyError:
        ch_num = 0
    
    if len(lonely_ch2) != ch2_num:
        return False
    elif len(lonely_ch) != ch_num:
        return False
    else:
        return True