"""Contribution model groups getter functions."""

from typing import Union

from rdkit import Chem

from ugropy.constants import (
    dort_ch2_hide,
    dort_ch_hide,
    dort_matrix,
    dort_problem,
    dort_subgroups,
    joback_ch2_hide,
    joback_ch_hide,
    joback_matrix,
    joback_problem,
    joback_subgroups,
    psrk_ch2_hide,
    psrk_ch_hide,
    psrk_matrix,
    psrk_problem,
    psrk_subgroups,
    unifac_ch2_hide,
    unifac_ch_hide,
    unifac_matrix,
    unifac_problem,
    unifac_subgroups,
)
from ugropy.core.get_model_groups import get_groups


def get_unifac_groups(
    identifier: Union[str, Chem.rdchem.Mol], identifier_type: str = "name"
) -> dict:
    """Get Classic LV-UNIFAC groups from molecule's name or SMILES.

    Parameters
    ----------
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
        Classic LV-UNIFAC subgroups.
    """

    unifac_groups = get_groups(
        unifac_subgroups,
        unifac_matrix,
        unifac_ch2_hide,
        unifac_ch_hide,
        unifac_problem,
        identifier,
        identifier_type,
    )

    return unifac_groups


def get_psrk_groups(
    identifier: Union[str, Chem.rdchem.Mol], identifier_type: str = "name"
) -> dict:
    """Get PSRK groups from molecule's name or SMILES.

    Parameters
    ----------
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
        PSRK subgroups.
    """
    psrk_groups = get_groups(
        psrk_subgroups,
        psrk_matrix,
        psrk_ch2_hide,
        psrk_ch_hide,
        psrk_problem,
        identifier,
        identifier_type,
    )

    return psrk_groups


def get_joback_groups(
    identifier: Union[str, Chem.rdchem.Mol], identifier_type: str = "name"
) -> dict:
    """Get Joback groups from molecule's name or SMILES.

    Parameters
    ----------
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
        Joback subgroups.
    """

    joback_groups = get_groups(
        joback_subgroups,
        joback_matrix,
        joback_ch2_hide,
        joback_ch_hide,
        joback_problem,
        identifier,
        identifier_type,
    )

    return joback_groups


def get_dortmund_groups(
    identifier: Union[str, Chem.rdchem.Mol], identifier_type: str = "name"
) -> dict:
    """Get Dortmund-UNIFAC groups from molecule's name or SMILES.

    Parameters
    ----------
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
        Dortmund-UNIFAC subgroups.
    """

    dortmund_groups = get_groups(
        dort_subgroups,
        dort_matrix,
        dort_ch2_hide,
        dort_ch_hide,
        dort_problem,
        identifier,
        identifier_type,
    )

    return dortmund_groups
