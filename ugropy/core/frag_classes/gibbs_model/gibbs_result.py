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
    calculate_r_q : bool
        If True, calculate R and Q values for the molecule.
    calculate_num_dict : bool, optional
        If True, calculate the subgroup numbers dictionary. Default is True.

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
    subgroups_num : dict
        Dictionary with the subgroup numbers and the number of times they
        appear in the molecule.
    """

    def __init__(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
        subgroups_info: pd.DataFrame,
        calculate_r_q: bool,
        calculate_num_dict: bool = True,
    ):
        super().__init__(molecule, subgroups, subgroups_atoms_indexes)

        r = 0.0
        q = 0.0

        # R and Q
        if self.subgroups != {}:
            if calculate_r_q:
                for group, n in self.subgroups.items():
                    r += n * subgroups_info.loc[group, "R"]
                    q += n * subgroups_info.loc[group, "Q"]

                self.r = r
                self.q = q
            else:
                self.r = None
                self.q = None
        else:
            self.r = None
            self.q = None

        # Subgroups numbers dictionary
        if self.subgroups != {}:
            if calculate_num_dict:
                self.subgroups_num = {}

                for group, occ in self.subgroups.items():
                    snum = int(subgroups_info.loc[group, "subgroup_number"])

                    self.subgroups_num[snum] = occ

            else:
                self.subgroups_num = None
