"""JobackFragmentationResult module."""

from typing import Union

import numpy as np

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors

from ugropy.constants import R
from ugropy.core.frag_classes.base.fragmentation_result import (
    FragmentationResult,
)


class JobackFragmentationResult(FragmentationResult):
    """Abdulelah Gani primary group contribution properties estimator.

    The class recieves either the Joback and Reid model's :cite:p:`joback1,
    joback2` groups, name or smiles of a molecule and estimates the its
    properties.

    Parameters
    ----------
    molecule : Chem.rdchem.Mol
        RDKit molecule object.
    subgroups : dict
        Dictionary of subgroups.
    subgroups_atoms_indexes : dict
        Dictionary of subgroups atoms indexes.
    info : pd.DataFrame
        Group's subgroups numbers.

    Attributes
    ----------
    subgroups : dict
        Abdulelah Gani primary functional groups of the molecule.
    subgroups_atoms_indexes : dict
        Abdulelah Gani primary functional groups atoms indexes.
    subgroups_numbers : dict
        Abdulelah Gani primary functional groups numbers.
    """

    def __init__(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
        info: pd.DataFrame,
    ) -> None:
        
