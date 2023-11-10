"""Groups module."""
from rdkit.Chem import Descriptors

from ugropy.core.model_getters import (
    get_psrk_groups,
    get_unifac_groups,
    instantiate_chem_object,
)


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
    psrk_groups : dict
        Predictive Soave-Redlich-Kwong subgroups.
    """

    def __init__(
        self,
        identifier: str,
        identifier_type: str = "name",
    ) -> None:
        self.identifier_type = identifier_type.lower()
        self.identifier = identifier
        self.chem_object = instantiate_chem_object(identifier, identifier_type)
        self.molecular_weight = Descriptors.MolWt(self.chem_object)

        # UNIFAC groups
        self.unifac_groups = get_unifac_groups(
            self.identifier, self.identifier_type
        )

        # PSRK groups
        self.psrk_groups = get_psrk_groups(
            self.identifier, self.identifier_type
        )

        # Joback
