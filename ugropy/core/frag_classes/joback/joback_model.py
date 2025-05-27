"""Joback fragmentation module."""

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
        the molecule).
    properties_contributions : pd.DataFrame, optional
        Group's properties contributions, by default None.

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Columns: 'smarts'
        (SMARTS representations of the group to detect its precense in the
        molecule).
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

        super().__init__(
            subgroups=subgroups,
            allow_overlapping=False,
            fragmentation_result=JobackFragmentationResult,
        )

        # Properties
        self.properties_contributions = properties_contributions

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
        search_nonoptimal: bool = False,
        solver_arguments: dict = {},
        normal_boiling_point: float = None,
    ) -> Union[JobackFragmentationResult, List[JobackFragmentationResult]]:
        """Get Jobacks groups from a molecule.

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
        normal_boiling_point : float, optional
            Experimental normal boiling point of the molecule on Kelvin. Its
            used to improve the properties calculations. Joback uses the
            estimated normal boiling point if no provided, by default None
        solver_arguments : dict, optional
            Dictionary with the arguments to be passed to the solver. For the
            DefaultSolver of ugropy you can change de PulP solver passing a
            dictionary like {"solver": "PULP_CBC_CMD"} and change the PulP
            solver. If empty it will use the default solver arguments, by
            default {}.

        Returns
        -------
        Union[JobackFragmentationResult, List[JobackFragmentationResult]]
            Fragmentation result. If search_multiple_solutions is False the
            return will be a FragmentationResult object, if True the return
            will be a list of FragmentationResult objects.
        """
        sol = super().get_groups(
            identifier=identifier,
            identifier_type=identifier_type,
            solver=solver,
            search_multiple_solutions=search_multiple_solutions,
            search_nonoptimal=search_nonoptimal,
            solver_arguments=solver_arguments,
            normal_boiling_point=normal_boiling_point,
            properties_contributions=self.properties_contributions,
        )

        return sol
