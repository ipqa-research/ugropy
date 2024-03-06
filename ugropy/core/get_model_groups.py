"""get_groups module.

Get the groups from a FragmentationModel.
"""

from typing import Union

import pandas as pd

from rdkit import Chem

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel

from .checks import (
    check_can_fit_atoms,
    check_has_composed,
    check_has_hiden,
    check_has_molecular_weight_right,
)
from .composed import correct_composed
from .detect_model_groups import detect_groups
from .fragmentation_object import Fragmentation
from .get_rdkit_object import instantiate_mol_object
from .problematics import correct_problematics


def get_groups(
    model: FragmentationModel,
    identifier: Union[str, Chem.rdchem.Mol],
    identifier_type: str = "name",
) -> Fragmentation:
    """Obtain the FragmentationModel's subgroups of an RDkit Mol object.

    Parameters
    ----------
    model: FragmentationModel
        FragmentationModel object.
    identifier : str or rdkit.Chem.rdchem.Mol
        Identifier of a molecule (name, SMILES or Chem.rdchem.Mol). Example:
        hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name, 'smiles' to provide the
        molecule SMILES representation or 'mol' to provide a
        rdkit.Chem.rdchem.Mol object, by default "name".

    Returns
    -------
    Fragmentation
        FragmentationModel's subgroups
    """
    # RDKit Mol object
    mol_object = instantiate_mol_object(identifier, identifier_type)

    # =========================================================================
    # Direct detection of groups presence and occurences
    # =========================================================================
    groups, groups_ocurrences = detect_groups(
        mol_object=mol_object, model=model
    )

    # =========================================================================
    # Filter the contribution matrix and sum over row to cancel the contribs
    # =========================================================================
    group_contributions = model.contribution_matrix.loc[groups, groups]
    group_contributions = group_contributions.mul(groups_ocurrences, axis=0)
    group_total_contributions = group_contributions.sum(axis=0)
    group_total_contributions.replace(0, pd.NA, inplace=True)

    mol_subgroups = group_total_contributions.dropna().to_dict()

    # =========================================================================
    # Check for the presence of problematic structures and correct.
    # =========================================================================
    mol_subgroups_corrected = correct_problematics(
        mol_object=mol_object,
        mol_subgroups=mol_subgroups,
        model=model,
    )

    # First exit
    if mol_subgroups_corrected == {}:
        # No functional groups were detected for the molecule. Example: H2O2
        return Fragmentation({}, mol_object, model)

    # =========================================================================
    # Check the presence of composed structures and check if the molecular
    # weight of the molecule is equals than the RDKit molecular weight.
    # =========================================================================
    right_mw = check_has_molecular_weight_right(
        mol_object=mol_object,
        mol_subgroups=mol_subgroups_corrected,
        model=model,
    )

    has_composed, _ = check_has_composed(
        mol_subgroups=mol_subgroups_corrected,
        model=model,
    )

    # =========================================================================
    # What to do according to the previous checks
    # =========================================================================
    if right_mw and not has_composed:
        # No need to do more, the solution was obtained.
        return Fragmentation(mol_subgroups_corrected, mol_object, model)
    elif not right_mw and not has_composed:
        # Nothing to do, the moelcule can't be modeled with FragmentationModel
        return Fragmentation({}, mol_object, model)
    elif not right_mw and has_composed:
        # Try fix the problem, the decomposition could still fail and return {}
        mol_subgroups_decomposed = correct_composed(
            mol_object=mol_object,
            mol_subgroups=mol_subgroups_corrected,
            model=model,
        )
        return Fragmentation(mol_subgroups_decomposed, mol_object, model)
    elif right_mw and has_composed:
        # Worst scenario, right mw and has composed, need check if has hidden
        has_hiden = check_has_hiden(mol_object, mol_subgroups_corrected, model)

        if has_hiden:
            mol_subgroups_decomposed = correct_composed(
                mol_object=mol_object,
                mol_subgroups=mol_subgroups_corrected,
                model=model,
            )

            fr = Fragmentation(mol_subgroups_decomposed, mol_object, model)

            return fr
        else:
            can_fit = check_can_fit_atoms(
                mol_object, mol_subgroups_corrected, model
            )

            if can_fit:
                fr = Fragmentation(mol_subgroups_corrected, mol_object, model)
                return fr
            else:
                mol_subgroups_decomposed = correct_composed(
                    mol_object=mol_object,
                    mol_subgroups=mol_subgroups_corrected,
                    model=model,
                )
                fr = Fragmentation(mol_subgroups_decomposed, mol_object, model)
                return fr
