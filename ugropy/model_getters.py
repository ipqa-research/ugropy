"""Contribution model groups getter functions."""
from typing import Union

import pubchempy as pcp

from rdkit import Chem

from ugropy.constants import (
    joback_ch2_hideouts,
    joback_ch_hideouts,
    joback_matrix,
    joback_problematics,
    joback_subgroups,
    problematic_structures,
    psrk_ch2_hideouts,
    psrk_ch_hideouts,
    psrk_matrix,
    psrk_subgroups,
    unifac_ch2_hideouts,
    unifac_ch_hideouts,
    unifac_matrix,
    unifac_subgroups,
)
from ugropy.core.get_groups import get_groups


def instantiate_chem_object(
    identifier: Union[str, Chem.rdchem.Mol], identifier_type: str = "name"
) -> Chem.rdchem.Mol:
    """Instantiate a rdkit Mol object from molecule's name or SMILES.

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
    rdkit.Chem.rdchem.Mol
        RDKit Mol object.
    """
    if identifier_type.lower() == "smiles":
        smiles = identifier
        chem_object = Chem.MolFromSmiles(smiles)

    elif identifier_type.lower() == "name":
        pcp_object = pcp.get_compounds(identifier, identifier_type)[0]
        smiles = pcp_object.canonical_smiles
        chem_object = Chem.MolFromSmiles(smiles)

    elif identifier_type.lower() == "mol":
        chem_object = identifier

        if not isinstance(chem_object, Chem.rdchem.Mol):
            raise ValueError(
                "If 'mol' identifier type is used, the identifier must be a "
                "rdkit.Chem.Chem.rdchem.Mol object."
            )

    else:
        raise ValueError(
            f"Identifier type: {identifier_type} not valid, use: 'name', "
            "'smiles' or 'mol'"
        )

    return chem_object


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
    chem_object = instantiate_chem_object(identifier, identifier_type)

    unifac_groups = get_groups(
        chem_object,
        unifac_subgroups,
        unifac_matrix,
        unifac_ch2_hideouts,
        unifac_ch_hideouts,
        problematic_structures,
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
    chem_object = instantiate_chem_object(identifier, identifier_type)

    psrk_groups = get_groups(
        chem_object,
        psrk_subgroups,
        psrk_matrix,
        psrk_ch2_hideouts,
        psrk_ch_hideouts,
        problematic_structures,
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
    chem_object = instantiate_chem_object(identifier, identifier_type)

    joback_groups = get_groups(
        chem_object,
        joback_subgroups,
        joback_matrix,
        joback_ch2_hideouts,
        joback_ch_hideouts,
        joback_problematics,
    )

    return joback_groups
