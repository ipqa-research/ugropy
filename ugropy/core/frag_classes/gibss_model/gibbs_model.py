"""GibbsModel fragmentation module."""

from collections import defaultdict
from typing import List, Union

import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.base.fragmentation_model import (
    FragmentationModel,
)
from ugropy.core.frag_classes.gibss_model.gibbs_result import (
    GibbsFragmentationResult,
)


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
        solutions_fragments: List[dict],
        search_multiple_solutions: bool,
    ) -> Union[GibbsFragmentationResult, List[GibbsFragmentationResult]]:
        """Get the solutions and return the GibbsFragmentationResult objects.

        Parameters
        ----------
        molecule : Chem.rdchem.Mol
            Rdkit mol object.
        solutions_fragments : List[dict]
            Fragments detected in the molecule.
        search_multiple_solutions : bool, optional
            Weather search for multiple solutions or not, by default False

        Returns
        -------
        Union[GibbsFragmentationResult, List[GibbsFragmentationResult]]
            List of GibbsFragmentationResult objects.
        """

        sols = []
        occurs = []

        for solution in solutions_fragments:
            occurrences = defaultdict(int)
            groups_atoms = defaultdict(list)

            for frag, atoms in solution.items():
                name = frag.split("_")[0]
                occurrences[name] += 1
                groups_atoms[name].append(atoms)

            if occurrences not in occurs:
                sols.append(
                    GibbsFragmentationResult(
                        molecule,
                        dict(occurrences),
                        dict(groups_atoms),
                        self.subgroups_info,
                    )
                )
                occurs.append(occurrences)

        if search_multiple_solutions:
            return sols
        else:
            return sols[0]
