"""get_groups module.

Get the groups from a FragmentationModel.
"""

from typing import Union

import pandas as pd

from rdkit import Chem

from .checks import (
    check_has_composed,
    check_has_hidden_ch2_ch,
    check_has_molecular_weight_right,
)
from .composed import correct_composed
from .detect_model_groups import detect_groups
from .get_rdkit_object import instantiate_mol_object
from .problematics import correct_problematics


def get_groups(
    subgroups: pd.DataFrame,
    subgroups_matrix: pd.DataFrame,
    ch2_hideouts: pd.DataFrame,
    ch_hideouts: pd.DataFrame,
    problematic_structures: pd.DataFrame,
    identifier: Union[str, Chem.rdchem.Mol],
    identifier_type: str = "name",
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
    identifier : str or rdkit.Chem.rdchem.Mol
        Identifier of a molecule (name, SMILES or Chem.rdchem.Mol). Example:
        hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name, 'smiles' to provide the
        molecule SMILES representation or 'mol' to provide a
        rdkit.Chem.rdchem.Mol object, by default "name".

    Returns
    -------
    dict
        FragmentationModel's subgroups
    """
    # RDKit Mol object
    mol_object = instantiate_mol_object(identifier, identifier_type)

    # Shorter names for DataFrames:
    df = subgroups.copy()
    dfm = subgroups_matrix.copy()
    dfp = problematic_structures.copy()

    # Groups found and the occurrence number of each one in chem_object.
    groups, groups_ocurrences = detect_groups(
        chem_object=mol_object, subgroups=df
    )

    # Filters the subgroups matrix into the detected groups only.
    dff = dfm.loc[groups][groups]

    # Multiply each group row by the occurrences of that group.
    dff = dff.mul(groups_ocurrences, axis=0)

    # Correction of problematic structures in dff.
    dff_corrected = correct_problematics(
        chem_object=mol_object,
        filtered_subgroups=dff,
        problematic_structures=dfp,
    )

    # import ipdb; ipdb.set_trace()

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
        chem_object=mol_object, chem_subgroups=chem_subgroups, subgroups=df
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
            chem_object=mol_object,
            chem_subgroups=chem_subgroups,
            subgroups=df,
            ch2_hideouts=ch2_hideouts,
            ch_hideouts=ch_hideouts,
        )
        return chem_subgroups
    elif right_mw and has_composed:
        has_hidden = check_has_hidden_ch2_ch(
            chem_object=mol_object,
            chem_subgroups=chem_subgroups,
            subgroups=subgroups,
            ch2_hideouts=ch2_hideouts,
            ch_hideouts=ch_hideouts,
        )
        if has_hidden:
            chem_subgroups = correct_composed(
                chem_object=mol_object,
                chem_subgroups=chem_subgroups,
                subgroups=df,
                ch2_hideouts=ch2_hideouts,
                ch_hideouts=ch_hideouts,
            )
            return chem_subgroups
        else:
            return chem_subgroups
