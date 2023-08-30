from rdkit import Chem
from rdkit.Chem import Descriptors

import numpy as np

import pandas as pd

from .detect_groups import group_matches
from ugropy.constants import ch2_hideouts, ch_hideouts


def check_has_molecular_weight_right(
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


def check_has_composed(
    chem_subgroups: dict,
    subgroups: pd.DataFrame,
) -> bool:
    composed_structures = subgroups[subgroups["composed"] == "y"].index
    
    for composed in composed_structures:
        if composed in chem_subgroups.keys():
            return True
        
    return False


def check_has_hidden_ch2_ch(
    chem_object: Chem.rdchem.Mol, 
    chem_subgroups: dict, 
    subgroups: pd.DataFrame
) -> bool:    
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
    