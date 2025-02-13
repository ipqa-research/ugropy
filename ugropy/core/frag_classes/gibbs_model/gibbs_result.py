"""Gibbs fragmentation result module."""

import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.base.fragmentation_result import (
    FragmentationResult,
)


class GibbsFragmentationResult(FragmentationResult):
    """Gibbs fragmentation result class.

    Parameters
    ----------
    molecule : Chem.rdchem.Mol
        RDKit molecule object.
    subgroups : dict
        Dictionary of subgroups.
    subgroups_atoms_indexes : dict
        Dictionary of subgroups atoms indexes.
    subgroups_info : pd.DataFrame
        DataFrame with subgroups information.

    Attributes
    ----------
    molecule : Chem.rdchem.Mol
        Molecule to fragment.
    subgroups : dict
        Dictionary with the subgroups and the number of times they appear in
        the molecule.
    subgroups_atoms : dict
        Dictionary with the subgroups and the atoms indexes that belong to
        each subgroup.
    r : float
        Gibss excess model R value estimation of the molecule.
    q : float
        Gibss excess model Q value estimation of the molecule.
    """

    def __init__(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
        subgroups_info: pd.DataFrame,
    ):
        super().__init__(molecule, subgroups, subgroups_atoms_indexes)

        r = 0.0
        q = 0.0

        if self.subgroups != {}:
            for group, n in self.subgroups.items():
                r += n * subgroups_info.loc[group, "R"]
                q += n * subgroups_info.loc[group, "Q"]

            self.r = r
            self.q = q
