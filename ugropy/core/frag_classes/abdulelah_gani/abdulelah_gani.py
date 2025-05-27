"""Abdulelah-Gani fragmentation module :cite:p:`gani`."""

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
    """Abdulelah-Gani fragmentation model :cite:p:`gani`.

    Parameters
    ----------
    abdulelah_gani_p : AbdulelahGaniPSTModel
        The primary Abdulelah-Gani fragmentation model.
    abdulelah_gani_s : AbdulelahGaniPSTModel
        The secondary Abdulelah-Gani fragmentation model.
    abdulelah_gani_t : AbdulelahGaniPSTModel
        The tertiary Abdulelah-Gani fragmentation model.
    properties_contributions : pd.DataFrame
        The contributions parameters of each group for each property of the
        model.
    properties_biases : pd.DataFrame
        The biases parameters of each property of the model.

    Attributes
    ----------
    primary_model : AbdulelahGaniPSTModel
        The primary Abdulelah-Gani fragmentation model.
    secondary_model : AbdulelahGaniPSTModel
        The secondary Abdulelah-Gani fragmentation model.
    tertiary_model : AbdulelahGaniPSTModel
        The tertiary Abdulelah-Gani fragmentation model.
    properties_contributions : pd.DataFrame
        The contributions parameters of each group for each property of the
        model.
    properties_biases : pd.DataFrame
        The biases parameters of each property of the model
    """

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
        search_nonoptimal: bool = False,
        solver_arguments: dict = {},
    ) -> Union[AGaniFragmentationResult, List[AGaniFragmentationResult]]:
        """Get the groups of the molecule.

        Parameters
        ----------
        identifier : Union[str, Chem.rdchem.Mol]
            Identifier of the molecule. You can use either the name of the
            molecule, the SMILEs of the molecule or a rdkit Mol object.
        identifier_type : str, optional
            Identifier type of the molecule. Use "name" if you are providing
            the molecules' name, "smiles" if you are providing the SMILES
            or "mol" if you are providing a rdkir mol object, by default "name"
        solver : ILPSolver, optional
            ILP solver class, by default DefaultSolver
        search_multiple_solutions : bool, optional
            Weather search for multiple solutions or not, by default False
            If False the return will be a FragmentationResult object, if True
            the return will be a list of FragmentationResult objects.
        search_nonoptimal : bool, optional
            If True, the solver will search for non-optimal solutions along
            with the optimal ones. This is useful when the user wants to find
            all possible combinations of fragments that cover the universe. By
            default False. If `search_multiple_solutions` is False, this
            parameter will be ignored.
        solver_arguments : dict, optional
            Dictionary with the arguments to be passed to the solver. For the
            DefaultSolver of ugropy you can change de PulP solver passing a
            dictionary like {"solver": "PULP_CBC_CMD"} and change the PulP
            solver. If empty it will use the default solver arguments, by
            default {}.

        Returns
        -------
        Union[AGaniFragmentationResult, List[AGaniFragmentationResult]]
            Fragmentation result. If search_multiple_solutions is False the
            return will be a FragmentationResult object, if True the return
            will be a list of FragmentationResult objects.
        """
        primary_groups = self.primary_model.get_groups(
            identifier,
            identifier_type,
            solver,
            search_multiple_solutions,
            search_nonoptimal,
            solver_arguments,
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
