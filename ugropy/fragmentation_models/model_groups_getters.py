"""Contribution model groups getter functions."""

from typing import Union

from rdkit import Chem

from ugropy.core.get_model_groups import get_groups
from ugropy.fragmentation_models.models import dortmund, joback, unifac, psrk


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

    unifac_groups = get_groups(unifac, identifier, identifier_type)

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
    psrk_groups = get_groups(psrk, identifier, identifier_type)

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

    joback_groups = get_groups(joback, identifier, identifier_type)

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

    dortmund_groups = get_groups(dortmund, identifier, identifier_type)

    return dortmund_groups
