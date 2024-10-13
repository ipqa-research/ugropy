"""FragmentationModel module.

All ugropy models (joback, unifac, psrk) are instances of the
FragmentationModule class.
"""

from collections import defaultdict

from typing import List, Union

import pandas as pd

from rdkit import Chem

import numpy as np

from ugropy.core.checks import check_atoms_fragments_presence
from ugropy.core.get_rdkit_object import instantiate_mol_object
from ugropy.core.frag_classes.base.fragmentation_result import (
    FragmentationResult,
)
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
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct).

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct).
    detection_mols : dict
        Dictionary cotaining all the rdkit Mol object from the detection_smarts
        subgroups column.
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        allow_overlapping: bool = False,
        check_molecular_weight: bool = False,
    ) -> None:
        self.subgroups = subgroups
        self.allow_overlapping = allow_overlapping
        self.check_molecular_weight = check_molecular_weight

        # Instantiate all de mol object from their smarts representation
        self.detection_mols = {}

        for group, row in self.subgroups.iterrows():
            self.detection_mols[group] = Chem.MolFromSmarts(row["smarts"])

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
    ) -> Union[FragmentationResult, List[FragmentationResult]]:
        # =====================================================================
        # Direct fragments detection
        # =====================================================================
        mol = instantiate_mol_object(identifier, identifier_type)

        detections = self.detect_fragments(mol)

        # No groups detected
        if detections == {}:
            return self.set_fragmentation_result(mol, [{}])

        # =====================================================================
        # Search for overlapping atoms and free atoms
        # =====================================================================
        overlapping_atoms, free_atoms = check_atoms_fragments_presence(
            mol, detections
        )

        # If there is free atoms in the molecule can't fragment with the model
        if np.size(free_atoms) > 0:
            return self.set_fragmentation_result(mol, [{}])

        # If no overlapping or the model allows overlapping we are done
        if np.size(overlapping_atoms) == 0 or self.allow_overlapping:
            return self.set_fragmentation_result(mol, [detections])

        # =====================================================================
        # Solve overlapping atoms
        # =====================================================================
        problem = solver(
            overlapping_atoms, detections, search_multiple_solutions
        )

        problem.solve()

        solutions = []

        for selection in problem.selected_fragments:
            solution = detections.copy()

            for frag in problem.overlapped_fragments:
                if frag not in selection:
                    solution.pop(frag)

            solutions.append(solution)

        return self.set_fragmentation_result(mol, solutions)

    def set_fragmentation_result(
        self,
        molecule: Chem.rdchem.Mol,
        solutions_fragments: List[dict],
    ) -> List[FragmentationResult]:

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
                    FragmentationResult(
                        molecule, dict(occurrences), dict(groups_atoms)
                    )
                )
                occurs.append(occurrences)

        return sols

    def detect_fragments(self, molecule: Chem.rdchem.Mol) -> dict:
        """Detect all the fragments in the molecule.

        Return a dictionary with the detected fragments as keys and a tuple
        with the atoms indexes of the fragment as values. For example, n-hexane
        for the UNIFAC model will return:

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
