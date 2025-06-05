"""FragmentationModel module.

All ugropy models (joback, unifac, psrk, etc) are instances of the
FragmentationModule class.
"""

from collections import defaultdict
from typing import List, Union

import numpy as np

import pandas as pd

from rdkit import Chem
from rdkit.Chem import rdMolEnumerator

from ugropy.core.checks import check_atoms_fragments_presence
from ugropy.core.frag_classes.base.fragmentation_result import (
    FragmentationResult,
)
from ugropy.core.get_rdkit_object import instantiate_mol_object
from ugropy.core.ilp_solvers.default_solver import DefaultSolver
from ugropy.core.ilp_solvers.ilp_solver import ILPSolver


class FragmentationModel:
    """FragmentationModel class.

    All ugropy supported models are an instance of this class. This class must
    be inherited to create a new type of FragmentationModel.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule).
    allow_overlapping : bool, optional
        Weather allow overlapping or not, by default False
    fragmentation_result : FragmentationResult, optional
        Fragmentation result class, by default FragmentationResult

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule).
    detection_mols : dict
        Dictionary cotaining all the rdkit Mol object from the detection_smarts
        subgroups column.
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        allow_overlapping: bool = False,
        allow_free_atoms: bool = False,
        fragmentation_result: FragmentationResult = FragmentationResult,
    ) -> None:
        self.subgroups = subgroups
        self.allow_overlapping = allow_overlapping
        self.allow_free_atoms = allow_free_atoms
        self.fragmentation_result = fragmentation_result

        # Instantiate all de mol object from their SMARTS representation
        self.detection_mols = {}

        for group, row in self.subgroups.iterrows():
            smarts = row["smarts"]

            if "LN:" in smarts:
                # Linknode expression from Chemaxon SMARTS
                smarts = smarts.replace("LN:", "|LN:") + "|"
                query = Chem.MolFromSmarts(smarts)
                enumerated_smarts = rdMolEnumerator.Enumerate(query)
                self.detection_mols[group] = enumerated_smarts
            else:
                self.detection_mols[group] = Chem.MolFromSmarts(smarts)

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
        search_nonoptimal: bool = False,
        solver_arguments: dict = {},
        **kwargs,
    ) -> Union[FragmentationResult, List[FragmentationResult]]:
        """Get the groups of a molecule.

        Parameters
        ----------
        identifier : Union[str, Chem.rdchem.Mol]
            Identifier of the molecule. You can use either the name of the
            molecule, the SMILEs of the molecule or a rdkit Mol object.
        identifier_type : str, optional
            Identifier type of the molecule. Use "name" if you are providing
            the molecules' name, "smiles" if you are providing the SMILES or
            "mol" if you are providing a rdkir mol object, by default "name"
        solver : ILPSolver, optional
            ILP solver class, by default DefaultSolver
        search_multiple_solutions : bool, optional
            Weather search for multiple solutions or not, by default False If
            False the return will be a FragmentationResult object, if True the
            return will be a list of FragmentationResult objects.
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
        Union[FragmentationResult, List[FragmentationResult]]
            Fragmentation result. If search_multiple_solutions is False the
            return will be a FragmentationResult object, if True the return
            will be a list of FragmentationResult objects.
        """
        # =====================================================================
        # Direct fragments detection
        # =====================================================================
        mol = self.mol_preprocess(
            instantiate_mol_object(identifier, identifier_type)
        )

        detections = self.detect_fragments(mol)

        # No groups detected
        if detections == {}:
            return self.set_fragmentation_result(
                mol, [{}], search_multiple_solutions, **kwargs
            )

        # =====================================================================
        # Search for overlapping atoms and free atoms
        # =====================================================================
        overlapping_atoms, free_atoms = check_atoms_fragments_presence(
            mol, detections
        )

        # If there is free atoms in the molecule can't fragment with the model
        if np.size(free_atoms) > 0 and not self.allow_free_atoms:
            return self.set_fragmentation_result(
                mol, [{}], search_multiple_solutions, **kwargs
            )

        # If no overlapping or the model allows overlapping we are done
        if np.size(overlapping_atoms) == 0 or self.allow_overlapping:
            return self.set_fragmentation_result(
                mol, [detections], search_multiple_solutions, **kwargs
            )

        # =====================================================================
        # Solve overlapping atoms
        # =====================================================================
        problem: ILPSolver = solver(
            overlapping_atoms,
            detections,
            search_multiple_solutions,
            search_nonoptimal,
            solver_arguments,
        )

        problem.solve()

        if not problem.selected_fragments:
            # This could happend, no solution found. Example:
            # "CC(C)(C)OC(=O)OC1=CC=CC=C1" on UNIFAC.
            return self.set_fragmentation_result(
                mol, [{}], search_multiple_solutions, **kwargs
            )

        solutions = []

        for selection in problem.selected_fragments:
            solution = detections.copy()

            for frag in problem.overlapped_fragments:
                if frag not in selection:
                    solution.pop(frag)

            solutions.append(solution)

        return self.set_fragmentation_result(
            mol, solutions, search_multiple_solutions, **kwargs
        )

    def mol_preprocess(self, mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
        """Preprocess the molecule.

        This method is called before the detection of the fragments. It can be
        used to preprocess the molecule before the detection of the fragments.
        This allow to use RDKit functions to preprocess the molecule and make
        your SMARTS detection easier. The default implementation does nothing
        and leave the mol object as it is. You can inherit this class and
        override this method to preprocess the molecule.

        Parameters
        ----------
        mol : Chem.rdchem.Mol
            Molecule to preprocess.

        Returns
        -------
        Chem.rdchem.Mol
            Preprocessed molecule.
        """
        return mol

    def set_fragmentation_result(
        self,
        molecule: Chem.rdchem.Mol,
        solutions_fragments: List[dict],
        search_multiple_solutions: bool = False,
        **kwargs,
    ) -> Union[FragmentationResult, List[FragmentationResult]]:
        """Process the solutions and return the FragmentationResult objects.

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
        Union[FragmentationResult, List[FragmentationResult]]
            List of FragmentationResult objects.
        """
        sols = []
        occurs = []  # solutions registry to avoid repeated solutions

        for solution in solutions_fragments:
            occurrences = defaultdict(int)
            groups_atoms = defaultdict(list)

            for frag, atoms in solution.items():
                name = frag.split("_")[0]
                occurrences[name] += 1
                groups_atoms[name].append(atoms)

            if occurrences not in occurs:
                sols.append(
                    self.fragmentation_result(
                        molecule,
                        dict(occurrences),
                        dict(groups_atoms),
                        **kwargs,
                    )
                )
                occurs.append(occurrences)

        if search_multiple_solutions:
            return sols
        else:
            return sols[0]

    def detect_fragments(self, molecule: Chem.rdchem.Mol) -> dict:
        """Detect all the fragments in the molecule.

        Return a dictionary with the detected fragments as keys and a tuple
        with the atoms indexes of the fragment as values. For example, n-hexane
        for the UNIFAC model will return:

        .. code-block:: python

            {
                'CH3_0': (0,),
                'CH3_1': (5,),
                'CH2_0': (1,),
                'CH2_1': (2,),
                'CH2_2': (3,),
                'CH2_3': (4,)
            }

        You may note that multiple occurrence of a fragment name will be
        indexed. The convention is always: <fragment_name>_i where `i` is the
        index of the occurrence.

        Parameters
        ----------
        mol : Chem.rdchem.Mol
            Molecule to detect the fragments.

        Returns
        -------
        dict
            Detected fragments in the molecule.
        """
        detected_fragments = {
            f"{fragment_name}_{i}": atoms_tuple
            for fragment_name, mol in self.detection_mols.items()
            for i, atoms_tuple in enumerate(molecule.GetSubstructMatches(mol))
        }

        return detected_fragments
