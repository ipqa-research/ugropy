from rdkit import Chem

import pandas as pd

from .detect_groups import detect_groups
from .correct_problematics import correct_problematics
from .checks import check_molecular_weight
from .correct_composed import correct_composed


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
    check_mw = check_molecular_weight(
        chem_object=chem_object, 
        chem_subgroups=chem_subgroups, 
        subgroups=df
    )
    
    if check_mw:
        return chem_subgroups
    else:
        chem_subgroups = correct_composed(
            chem_object=chem_object,
            chem_subgroups=chem_subgroups,
            subgroups=df
        )
        return chem_subgroups


