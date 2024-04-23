"""get_rdkit_object module."""

from functools import cache
from typing import Union

import pubchempy as pcp

from rdkit import Chem


@cache
def instantiate_mol_object(
    identifier: Union[str, Chem.rdchem.Mol], identifier_type: str = "name"
) -> Chem.rdchem.Mol:
    """Instantiate a RDKit Mol object from molecule's name or SMILES.

    Parameters
    ----------
    identifier : str or rdkit.Chem.rdchem.Mol
        Identifier of a molecule (name, SMILES or rdkit.Chem.rdchem.Mol).
        Example: hexane or CCCCCC for name or SMILES respectively.
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
