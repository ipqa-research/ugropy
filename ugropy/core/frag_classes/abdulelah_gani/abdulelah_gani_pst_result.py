"""AganiFragmentationResult module."""

import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.base.fragmentation_result import (
    FragmentationResult,
)


class AGaniPSTFragmentationResult(FragmentationResult):
    """Abdulelah-Gani group contribution properties estimator.

    This class is to instantiate the Abdulelah-Gani primary, secondary and
    tertiary models.

    Parameters
    ----------
    molecule : Chem.rdchem.Mol
        RDKit molecule object.
    subgroups : dict
        Dictionary of subgroups.
    subgroups_atoms_indexes : dict
        Dictionary of subgroups atoms indexes.
    subgroups_info : pd.DataFrame
        Group's subgroups numbers.

    Attributes
    ----------
    subgroups : dict
        Abdulelah-Gani groups of the molecule.
    subgroups_atoms_indexes : dict
        Abdulelah-Gani groups atoms indexes.
    subgroups_numbers : dict
        Abdulelah-Gani groups numbers.
    """

    def __init__(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
        subgroups_info: pd.DataFrame,
    ) -> None:

        super().__init__(molecule, subgroups, subgroups_atoms_indexes)

        self.subgroups_numbers = {}

        if self.subgroups != {}:
            for group, n in self.subgroups.items():
                group_number = int(subgroups_info.loc[group, "group_number"])
                self.subgroups_numbers[group_number] = n
