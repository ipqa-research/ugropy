"""get_rdkit_object module."""

from functools import cache
from typing import Union

import pubchempy as pcp

from rdkit import Chem


def instantiate_mol_object(
    identifier: Union[str, Chem.rdchem.Mol], identifier_type: str = "name"
) -> Chem.rdchem.Mol:
    """Instantiate a RDKit Mol object from molecule's name, SMILES or mol.

    In case that the input its already a RDKit Mol object, the function will
    return the input object.

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
        chem_object = instantiate_mol_from_name(identifier)

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


@cache
def instantiate_mol_from_name(name: str) -> Chem.rdchem.Mol:
    """Instantiate a RDKit Mol object from molecule's name.

    The funcion uses `pubchempy` to get the molecule's SMILES representation
    from the molecule's name and then instantiate a RDKit Mol object.

    Parameters
    ----------
    name : str
        Name of the molecule.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        RDKit Mol object.
    """
    try:
        pcp_object = pcp.get_compounds(name, "name")[0]
        smiles = pcp_object.canonical_smiles
        chem_object = Chem.MolFromSmiles(smiles)
    except IndexError:
        raise ValueError(
            f"Could not find a molecule with the name '{name}' on " "PubChem"
        )

    return chem_object
