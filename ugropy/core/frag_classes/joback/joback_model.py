"""GibbsModel fragmentation module."""

from collections import defaultdict
from typing import List, Union

import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.base.fragmentation_model import (
    FragmentationModel,
)
from ugropy.core.frag_classes.joback.joback_result import (
    JobackFragmentationResult,
)
from ugropy.core.ilp_solvers.default_solver import DefaultSolver
from ugropy.core.ilp_solvers.ilp_solver import ILPSolver


class JobackModel(FragmentationModel):
    """Joback Fragmentation model dedicated to properties estimation models.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct).
    properties_contributions : pd.DataFrame, optional
        Group's properties contributions, by default None

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct).
    detection_mols : dict
        Dictionary cotaining all the rdkit Mol object from the detection_smarts
        subgroups column
    properties_contributions : pd.DataFrame
        Group's properties contributions.
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        properties_contributions: Union[pd.DataFrame, None],
    ) -> None:

        super().__init__(subgroups)

        # Properties
        self.properties_contributions = properties_contributions

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
        normal_boiling_point: float = None,
    ) -> Union[JobackFragmentationResult, List[JobackFragmentationResult]]:

        sol = super().get_groups(
            identifier=identifier,
            identifier_type=identifier_type,
            solver=solver,
            search_multiple_solutions=search_multiple_solutions,
            normal_boiling_point=normal_boiling_point,
        )

        return sol

    def set_fragmentation_result(
        self,
        molecule: Chem.rdchem.Mol,
        solutions_fragments: List[dict],
        normal_boiling_point: float,
    ) -> List[JobackFragmentationResult]:

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
                    JobackFragmentationResult(
                        molecule,
                        dict(occurrences),
                        dict(groups_atoms),
                        self.properties_contributions,
                        normal_boiling_point=normal_boiling_point,
                    )
                )
                occurs.append(occurrences)

        return sols
