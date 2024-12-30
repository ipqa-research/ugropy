"""Joback fragmentation module."""

from typing import List, Union

import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_pst import (
    AbdulelahGaniPSTModel,
)
from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_result import (
    AGaniFragmentationResult,
)
from ugropy.core.ilp_solvers.default_solver import DefaultSolver
from ugropy.core.ilp_solvers.ilp_solver import ILPSolver


class AbdulelahGaniModel:
    def __init__(
        self,
        abdulelah_gani_p: AbdulelahGaniPSTModel,
        abdulelah_gani_s: AbdulelahGaniPSTModel,
        abdulelah_gani_t: AbdulelahGaniPSTModel,
        properties_contributions: pd.DataFrame,
        properties_biases: pd.DataFrame,
    ) -> None:

        self.primary_model = abdulelah_gani_p
        self.secondary_model = abdulelah_gani_s
        self.tertiary_model = abdulelah_gani_t
        self.properties_contributions = properties_contributions
        self.properties_biases = properties_biases

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
    ) -> Union[AGaniFragmentationResult, List[AGaniFragmentationResult]]:

        primary_groups = self.primary_model.get_groups(
            identifier, identifier_type, solver, search_multiple_solutions
        )

        secondary_groups = self.secondary_model.get_groups(
            identifier,
            identifier_type,
            solver,
            search_multiple_solutions=False,
        )

        tertiary_groups = self.tertiary_model.get_groups(
            identifier,
            identifier_type,
            solver,
            search_multiple_solutions=False,
        )

        if not search_multiple_solutions:
            result = AGaniFragmentationResult(
                primary_groups.molecule,
                primary_groups,
                secondary_groups,
                tertiary_groups,
                self.properties_contributions,
                self.properties_biases,
            )
        else:
            result = []

            for psol in primary_groups:
                result.append(
                    AGaniFragmentationResult(
                        psol.molecule,
                        psol,
                        secondary_groups,
                        tertiary_groups,
                        self.properties_contributions,
                        self.properties_biases,
                    )
                )

        return result
