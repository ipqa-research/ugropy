"""Joback fragmentation module."""

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
        normal_boiling_point : float, optional
            Experimental normal boiling point of the molecule on Kelvin. Its
            used to improve the properties calculations. Joback uses the
            estimated normal boiling point if no provided, by default None

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
            normal_boiling_point=normal_boiling_point,
        )

        return sol

    def set_fragmentation_result(
        self,
        molecule: Chem.rdchem.Mol,
        solutions_fragments: List[dict],
        search_multiple_solutions: bool,
        normal_boiling_point: float,
    ) -> Union[JobackFragmentationResult, List[JobackFragmentationResult]]:
        """Get the solutions and return the JobackFragmentationResult objects.

        Parameters
        ----------
        molecule : Chem.rdchem.Mol
            Rdkit mol object.
        solutions_fragments : List[dict]
            Fragments detected in the molecule.
        search_multiple_solutions : bool, optional
            Weather search for multiple solutions or not, by default False
        normal_boiling_point : float
            Experimental normal boiling point of the molecule on Kelvin.

        Returns
        -------
        Union[JobackFragmentationResult, List[JobackFragmentationResult]]
            Fragmentation result. If search_multiple_solutions is False the
            return will be a FragmentationResult object, if True the return
            will be a list of FragmentationResult objects.
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
                    JobackFragmentationResult(
                        molecule,
                        dict(occurrences),
                        dict(groups_atoms),
                        self.properties_contributions,
                        normal_boiling_point=normal_boiling_point,
                    )
                )
                occurs.append(occurrences)

        if search_multiple_solutions:
            return sols
        else:
            return sols[0]
