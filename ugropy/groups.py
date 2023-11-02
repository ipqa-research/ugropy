"""Groups module."""
from ugropy.core.model_getters import (
    instantiate_chem_object,
    get_psrk_groups,
    get_unifac_groups,
)
from rdkit.Chem import Descriptors


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
    dortmund : bool, optional
        If True the algorithm will try to get the Dortmund groups. If False 
        this will be skiped, by default "True".

    Attributes
    ----------
    identifier : str
        Identifier of a molecule. Example: hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name or 'smiles' to provide the
        molecule SMILES representation, by default "name".
    chem_object : rdkit.Chem.rdchem.Mol
        RDKit Mol object.
    molecular_weight : float
        Molecule's molecular weight from rdkit.Chem.Descriptors.MolWt [g/mol].
    unifac_groups : dict
        Classic LV-UNIFAC subgroups.
    """

    def __init__(
        self,
        identifier: str,
        identifier_type: str = "name",
        unifac: bool = True,
        psrk: bool = True,
        dortmund: bool = True,
    ) -> None:
        self.identifier_type = identifier_type.lower()
        self.identifier = identifier
        self.chem_object = instantiate_chem_object(identifier, identifier_type)
        self.molecular_weight = Descriptors.MolWt(self.chem_object)

        # =====================================================================
        # UNIFAC groups
        # =====================================================================
        if unifac:
            self.unifac_groups = get_unifac_groups(
                self.identifier, self.identifier_type
            )
        else:
            self.unifac_groups = {}

        # =====================================================================
        # PSRK groups
        # =====================================================================
        if psrk:
            self.psrk_groups = get_psrk_groups(
                self.identifier, self.identifier_type
            )
        else:
            self.psrk_groups = {}

        # =====================================================================
        # Dortmund groups
        # =====================================================================
        if dortmund:
            self.dortmund_groups = get_groups(
                self.chem_object,
                dortmund_subgroups,
                dortmund_matrix,
                psrk_ch2_hideouts,
                psrk_ch_hideouts,
                problematic_structures,
            )
        else:
            self.dortmund_groups = {}