"""FragmentationModel module.

All ugropy models (joback, unifac, psrk) are instances of the
FragmentationModule class.
"""
from abc import ABC, abstractmethod

import pandas as pd

from rdkit import Chem


class FragmentationModel(ABC):
    """FragmentationModel class.

    All ugropy supported models are an instance of this class. This class must
    be inherited to create a new type of FragmentationModel.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct).

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct).
    detection_mols : dict
        Dictionary cotaining all the rdkit Mol object from the detection_smarts
        subgroups column.
    """

    def __init__(self, subgroups: pd.DataFrame) -> None:
        self.subgroups = subgroups
        
        # Instantiate all de mol object from their smarts representation
        detection_mols = {}

        for group, row in self.subgroups.iterrows():
            detection_mols[group] = Chem.MolFromSmarts(row["smarts"])

    
    def detect_groups(self, molecule: Chem.Mol) -> pd.DataFrame:
        """Detect all the groups in the molecule.

        Return a dictionary with the detected groups as keys and a tuple of
        tuples containing the molecule's atoms that participate in the group
        occurrences.

        Parameters
        ----------
        mol : Chem.Mol
            Molecule to detect the groups.

        Returns
        -------
        dict
            Detected groups in the molecule.
        """
        detected_groups = {}
        for group, mol in self.detection_mols.items():
            matches = molecule.GetSubstructMatches(mol)
                
            if matches:
                detected_groups[group] = matches

        return detected_groups
    
    @abstractmethod
    def set_fragmentation_result(
        self,
        molecule: Chem.Mol,
        subgroups_occurrences: dict,
        subgroups_atoms_indexes: dict,
    ) -> "FragmentationResult":
        
        raise NotImplementedError("Abstract Method not implemented.")
