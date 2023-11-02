"""Contribution model groups getter functions."""
import pubchempy as pcp

from rdkit import Chem

from ugropy.constants import (
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
    identifier: str, identifier_type: str
) -> Chem.rdchem.Mol:
    """Instantiates a rdkit Mol object from molecule's name or SMILES.

    Parameters
    ----------
    identifier : str
        Identifier of a molecule. Example: hexane or CCCCCC.
    identifier_type : str
        Use 'name' to search a molecule by name or 'smiles' to provide the
        molecule SMILES representation.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        RDKit Mol object.
    """
    if identifier_type.lower() == "smiles":
        smiles = identifier
    elif identifier_type.lower() == "name":
        pcp_object = pcp.get_compounds(identifier, identifier_type)[0]
        smiles = pcp_object.canonical_smiles
    else:
        raise ValueError

    chem_object = Chem.MolFromSmiles(smiles)

    return chem_object


def get_unifac_groups(identifier: str, identifier_type: str = "name") -> dict:
    """Get Classic LV-UNIFAC groups from molecule's name or SMILES.

    Parameters
    ----------
    identifier : str
        Identifier of a molecule. Example: hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name or 'smiles' to provide the
        molecule SMILES representation, by default "name".

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


def get_psrk_groups(identifier: str, identifier_type: str = "name") -> dict:
    """Get PSRK groups from molecule's name or SMILES.

    Parameters
    ----------
    identifier : str
        Identifier of a molecule. Example: hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name or 'smiles' to provide the
        molecule SMILES representation, by default "name".

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