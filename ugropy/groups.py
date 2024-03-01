"""Groups module."""

from rdkit.Chem import Descriptors

from ugropy.core.get_model_groups import get_groups
from ugropy.core.get_rdkit_object import instantiate_mol_object
from ugropy.fragmentation_models.models import psrk, unifac
from ugropy.properties.joback_properties import JobackProperties


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
    mol_object : rdkit.Chem.rdchem.Mol
        RDKit Mol object.
    molecular_weight : float
        Molecule's molecular weight from rdkit.Chem.Descriptors.MolWt [g/mol].
    unifac : Fragmentation
        Classic LV-UNIFAC subgroups.
    psrk : Fragmentation
        Predictive Soave-Redlich-Kwong subgroups.
    joback : JobackProperties
        JobackProperties object that contains the Joback subgroups and the
        estimated properties of the molecule.
    """

    def __init__(
        self,
        identifier: str,
        identifier_type: str = "name",
        normal_boiling_temperature: float = None,
    ) -> None:
        self.identifier_type = identifier_type.lower()
        self.identifier = identifier
        self.mol_object = instantiate_mol_object(identifier, identifier_type)
        self.molecular_weight = Descriptors.MolWt(self.mol_object)

        # UNIFAC
        self.unifac = get_groups(unifac, self.identifier, self.identifier_type)

        # PSRK
        self.psrk = get_groups(psrk, self.identifier, self.identifier_type)

        # Joback
        self.joback = JobackProperties(
            self.identifier, self.identifier_type, normal_boiling_temperature
        )
