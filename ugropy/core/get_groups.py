"""get_groups module.

Gets the UNIFAC's subgroups of an RDKit Mol object.
"""
import pandas as pd

from rdkit import Chem

from .checks import (
    check_has_composed,
    check_has_hidden_ch2_ch,
    check_has_molecular_weight_right,
)
from .correct_composed import correct_composed
from .correct_problematics import correct_problematics
from .detect_groups import detect_groups


def get_groups(
    chem_object: Chem.rdchem.Mol,
    subgroups: pd.DataFrame,
    subgroups_matrix: pd.DataFrame,
    ch2_hideouts: pd.DataFrame,
    ch_hideouts: pd.DataFrame,
    problematic_structures: pd.DataFrame,
):
    """Obtain the UNIFAC's model subgroups of chem_object.

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Mol object.
    subgroups : pd.DataFrame
        DataFrame that contains the subgroups, their SMARTS representation and
        subgroup count contribution.
    subgroups_matrix : pd.DataFrame
        DataFrame that cotains the subgroups contribution matrix.
    problematic_structures : pd.DataFrame
        DataFrame that contains the SMARTS representation of the problematic
        structures and their contribution correction.

    Returns
    -------
    dict
        UNIFAC's subgroups
    """
    # Shorter names for DataFrames:
    df = subgroups.copy()
    dfm = subgroups_matrix.copy()
    dfp = problematic_structures.copy()

    # Groups found and the occurrence number of each one in chem_object.
    groups, groups_ocurrences = detect_groups(
        chem_object=chem_object, subgroups=df
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
        # No functional groups were detected for the molecule. Example:
        # hydrogen peroxide.
        return chem_subgroups

    # Check for composed structures.
    right_mw = check_has_molecular_weight_right(
        chem_object=chem_object, chem_subgroups=chem_subgroups, subgroups=df
    )

    has_composed = check_has_composed(
        chem_subgroups=chem_subgroups, subgroups=subgroups
    )

    if right_mw and not has_composed:
        return chem_subgroups
    elif not right_mw and not has_composed:
        return {}
    elif not right_mw and has_composed:
        chem_subgroups = correct_composed(
            chem_object=chem_object,
            chem_subgroups=chem_subgroups,
            subgroups=df,
            ch2_hideouts=ch2_hideouts,
            ch_hideouts=ch_hideouts,
        )
        return chem_subgroups
    elif right_mw and has_composed:
        has_hidden = check_has_hidden_ch2_ch(
            chem_object=chem_object,
            chem_subgroups=chem_subgroups,
            subgroups=subgroups,
            ch2_hideouts=ch2_hideouts,
            ch_hideouts=ch_hideouts,
        )
        if has_hidden:
            chem_subgroups = correct_composed(
                chem_object=chem_object,
                chem_subgroups=chem_subgroups,
                subgroups=df,
                ch2_hideouts=ch2_hideouts,
                ch_hideouts=ch_hideouts,
            )
            return chem_subgroups
        else:
            return chem_subgroups
