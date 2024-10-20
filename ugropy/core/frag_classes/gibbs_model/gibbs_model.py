"""GibbsModel fragmentation module."""

from collections import defaultdict
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

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
        **kwargs,
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
            **kwargs,
        )

        return sol

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
