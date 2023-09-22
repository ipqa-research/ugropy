"""Groups module."""
import pubchempy as pcp

from rdkit import Chem

from .constants import (
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
        molecule SMILES representation, by default "name".
    unifac : bool, optional
        If True the algorithm will try to get the classic LV-UNIFAC groups. If
        False this will be skiped, by default "True".
    psrk : bool, optional
        If True the algorithm will try to get the PSRK groups. If False this
        will be skiped, by default "True".

    Attributes
    ----------
    identifier : str
        Identifier of a molecule. Example: hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name or 'smiles' to provide the
        molecule SMILES representation, by default "name".
    chem_object : rdkit.Chem.rdchem.Mol
        RDKit Mol object.
    unifac_groups : dict
        Classic LV-UNIFAC subgroups.
    """

    def __init__(
        self,
        identifier: str,
        identifier_type: str = "name",
        unifac: bool = True,
        psrk: bool = True,
    ) -> None:
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

        if unifac:
            self.unifac_groups = get_groups(
                self.chem_object,
                unifac_subgroups,
                unifac_matrix,
                unifac_ch2_hideouts,
                unifac_ch_hideouts,
                problematic_structures,
            )
        else:
            self.unifac_groups = {}

        if psrk:
            self.psrk_groups = get_groups(
                self.chem_object,
                psrk_subgroups,
                psrk_matrix,
                psrk_ch2_hideouts,
                psrk_ch_hideouts,
                problematic_structures,
            )
        else:
            self.psrk_groups = {}
