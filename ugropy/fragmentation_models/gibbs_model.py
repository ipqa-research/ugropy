"""GibbsModel fragmentation module."""

from typing import Union

import pandas as pd

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel
from ugropy.fragmentation_models.frag_result import FragmentationResult

from rdkit import Chem


class GibbsModel(FragmentationModel):
    """GibbsModel it's a fragmentation model dedicated to Gibbs excess models.

    unifac, psrk, dortmund are instances of this class.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct)
    subgroups_info : Union[pd.DataFrame, None], optional
        Information of the model's subgroups (R, Q, subgroup_number,
        main_group), by default None

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct)
    detection_mols : dict
        Dictionary cotaining all the rdkit Mol object from the detection_smarts
        subgroups column
    subgroups_info : pd.DataFrame
        Information of the model's subgroups. Columns: R, Q, subgroup_number,
        main_group. Index: 'group' (subgroups names)
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        subgroups_info: Union[pd.DataFrame, None] = None,
    ) -> None:
        super().__init__(subgroups)

        # subgroups info
        if subgroups_info is None:
            self.subgroups_info = pd.DataFrame(
                [],
                columns=["group", "subgroup_number", "main_group", "R", "Q"],
            ).set_index("group")
        else:
            self.subgroups_info = subgroups_info

    def set_fragmentation_result(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
    ) -> "GibbsFragmentationResult":

        result = GibbsFragmentationResult(
            self, molecule, subgroups, subgroups_atoms_indexes
        )

        return result


class GibbsFragmentationResult(FragmentationResult):
    def __init__(
        self,
        model: GibbsModel,
        molecule: Chem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
    ):
        super().__init__(molecule, subgroups, subgroups_atoms_indexes)

        self.r = 0.0
        self.q = 0.0

        self._set_r_and_q(model)

    def _set_r_and_q(self, model: GibbsModel) -> None:
        """Set R and Q values to the subgroups."""
        r = 0.0
        q = 0.0

        for group, n in self.subgroups.items():
            r += n * model.subgroups_info.loc[group, "R"]
            q += n * model.subgroups_info.loc[group, "Q"]

        self.r = r
        self.q = q
