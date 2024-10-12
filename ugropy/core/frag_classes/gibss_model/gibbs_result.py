from rdkit import Chem

from ugropy.core.frag_classes.gibss_model.gibbs_model import GibbsModel
from ugropy.core.frag_classes.base.fragmentation_result import (
    FragmentationResult,
)


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
