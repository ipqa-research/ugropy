"""Groups module."""
import pubchempy as pcp

from rdkit import Chem

from .constants import problematic_structures, unifac_matrix, unifac_subgroups
from .core.get_groups import get_groups


class Groups:
    """Group class.

    Stores the solved UNIFAC subgroups of a molecule.

    Parameters
    ----------
    identifier : str
        Identifier of a molecule. Example: hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name or 'smiles' to provide the
        molecule SMILES representation, by default "name"

    Attributes
    ----------
    identifier : str
        Identifier of a molecule. Example: hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name or 'smiles' to provide the
        molecule SMILES representation, by default "name"
    chem_object : rdkit.Chem.rdchem.Mol
        RDKit Mol object.
    unifac_groups : dict
        Classic LV-UNIFAC subgroups.
    """

    def __init__(self, identifier: str, identifier_type: str = "name") -> None:
        self.identifier = identifier.lower()
        self.identifier_type = identifier_type.lower()

        if self.identifier_type == "smiles":
            self.smiles = identifier
            self.chem_object = Chem.MolFromSmiles(self.smiles)
        else:
            pcp_object = pcp.get_compounds(
                self.identifier, self.identifier_type
            )[0]
            self.smiles = pcp_object.canonical_smiles
            self.chem_object = Chem.MolFromSmiles(self.smiles)

        self.unifac_groups = get_groups(
            self.chem_object,
            unifac_subgroups,
            unifac_matrix,
            problematic_structures,
        )
