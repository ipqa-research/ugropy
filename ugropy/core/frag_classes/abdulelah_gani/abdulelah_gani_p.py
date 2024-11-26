"""Joback fragmentation module."""

from collections import defaultdict
from typing import List, Union

import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_p_result import (
    AGaniFragmentationResult,
)
from ugropy.core.frag_classes.base.fragmentation_model import (
    FragmentationModel,
)


class AbdulelahGaniPrimaryModel(FragmentationModel):
    """Abdulelah-Gani model dedicated to properties estimation models.

    Class to construct the primary structures detector for the Abdulelah-Gani
    properties estimation model :cite:p:`gani`.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule).
    info : pd.DataFrame
        Group's subgroups numbers.

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule).
    detection_mols : dict
        Dictionary cotaining all the rdkit Mol object from the detection_smarts
        subgroups column
    info : pd.DataFrame
        Group's subgroups numbers.
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        info: pd.DataFrame,
    ) -> None:

        super().__init__(subgroups)
        self.info = info

    def set_fragmentation_result(
        self,
        molecule: Chem.rdchem.Mol,
        solutions_fragments: List[dict],
        search_multiple_solutions: bool,
    ) -> Union[AGaniFragmentationResult, List[AGaniFragmentationResult]]:
        """Get the solutions and return the AGaniFragmentationResult objects.

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
        Union[AGaniFragmentationResult, List[AGaniFragmentationResult]]
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
                    AGaniFragmentationResult(
                        molecule,
                        dict(occurrences),
                        dict(groups_atoms),
                        self.info,
                    )
                )
                occurs.append(occurrences)

        if search_multiple_solutions:
            return sols
        else:
            return sols[0]
