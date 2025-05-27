"""GibbsModel fragmentation module."""

from typing import List, Union

import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.base.fragmentation_model import (
    FragmentationModel,
)
from ugropy.core.frag_classes.gibbs_model.gibbs_result import (
    GibbsFragmentationResult,
)
from ugropy.core.ilp_solvers.default_solver import DefaultSolver
from ugropy.core.ilp_solvers.ilp_solver import ILPSolver


class GibbsModel(FragmentationModel):
    """GibbsModel it's a fragmentation model dedicated to Gibbs excess models.

    unifac, psrk, dortmund are instances of this class.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule).
    subgroups_info : Union[pd.DataFrame, None], optional
        Information of the model's subgroups (R, Q, subgroup_number,
        main_group), by default None

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Columns: 'smarts'
        (SMARTS representations of the group to detect its precense in the
        molecule).
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
        calculate_r_q: bool = True,
    ) -> None:
        super().__init__(
            subgroups=subgroups,
            allow_overlapping=False,
            fragmentation_result=GibbsFragmentationResult,
        )

        self._calculate_r_q = calculate_r_q

        # subgroups info
        if subgroups_info is None:
            self.subgroups_info = pd.DataFrame(
                [],
                columns=["group", "subgroup_number", "main_group", "R", "Q"],
            ).set_index("group")
        else:
            self.subgroups_info = subgroups_info

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
        search_nonoptimal: bool = False,
        solver_arguments: dict = {},
    ) -> Union[GibbsFragmentationResult, List[GibbsFragmentationResult]]:
        """Get the groups of a molecule.

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
        Union[GibbsFragmentationResult, List[GibbsFragmentationResult]]
            Fragmentation result. If search_multiple_solutions is False the
            return will be a FragmentationResult object, if True the return
            will be a list of FragmentationResult objects.
        """
        sol = super().get_groups(
            identifier,
            identifier_type,
            solver,
            search_multiple_solutions,
            search_nonoptimal,
            solver_arguments,
            subgroups_info=self.subgroups_info,
            calculate_r_q=self._calculate_r_q,
        )

        return sol
