from typing import List, Union

from rdkit import Chem
from rdkit.Chem import Descriptors

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


class Fragmentation:
    def __init__(
        self,
        mol_subgroups: Union[dict, List[dict]],
        mol_object: Chem.rdchem.Mol,
        model: FragmentationModel,
    ) -> None:
        self.subgroups = mol_subgroups
        self.mol_object = mol_object
        self.molecular_weight = Descriptors.MolWt(self.mol_object)
        self.smiles = Chem.MolToSmiles(self.mol_object)
        self._model = model
