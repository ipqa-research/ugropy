import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.base.fragmentation_result import (
    FragmentationResult,
)


class GibbsFragmentationResult(FragmentationResult):
    def __init__(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
        subgroups_info: pd.DataFrame,
    ):
        super().__init__(molecule, subgroups, subgroups_atoms_indexes)

        r = 0.0
        q = 0.0

        if self.subgroups != {}:
            for group, n in self.subgroups.items():
                r += n * subgroups_info.loc[group, "R"]
                q += n * subgroups_info.loc[group, "Q"]

            self.r = r
            self.q = q
