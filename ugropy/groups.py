"""Groups module."""

from rdkit.Chem import Descriptors

from ugropy.core.get_model_groups import get_groups
from ugropy.core.get_rdkit_object import instantiate_mol_object
from ugropy.fragmentation_models.models import psrk, unifac
from ugropy.joback_properties import Joback


class Groups:
    """Group class.

    Stores the solved FragmentationModels subgroups of a molecule.

    Parameters
    ----------
    identifier : str or rdkit.Chem.rdchem.Mol
        Identifier of a molecule (name, SMILES or Chem.rdchem.Mol). Example:
        hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name, 'smiles' to provide the
        molecule SMILES representation or 'mol' to provide a
        rdkit.Chem.rdchem.Mol object, by default "name".
    normal_boiling_temperature : float, optional
        If provided, will be used to estimate critical temperature, acentric
        factor, and vapor pressure instead of the estimated normal boiling
        point in the Joback group contribution model, by default None.

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
    joback : Joback
        Joback object that contains the Joback subgroups and the estimated
        properties of the molecule.
    """

    def __init__(
        self,
        identifier: str,
        identifier_type: str = "name",
        normal_boiling_temperature: float = None,
    ) -> None:
        self.identifier_type = identifier_type.lower()
        self.identifier = identifier
        self.chem_object = instantiate_mol_object(identifier, identifier_type)
        self.molecular_weight = Descriptors.MolWt(self.chem_object)

        # UNIFAC groups
        self.unifac_groups = get_groups(
            unifac, self.identifier, self.identifier_type
        )

        # PSRK groups
        self.psrk_groups = get_groups(
            psrk, self.identifier, self.identifier_type
        )

        # Joback
        self.joback = Joback(
            self.identifier, self.identifier_type, normal_boiling_temperature
        )
