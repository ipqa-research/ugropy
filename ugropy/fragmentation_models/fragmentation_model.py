"""FragmentationModel module.

All ugropy models (joback, unifac, psrk) are instances of the
FragmentationModule class.
"""

import json
from typing import List, Union

import numpy as np

import pandas as pd

from rdkit import Chem


class FragmentationModel:
    """FragmentationModel class.

    All ugropy supported models are an instance of this class. This class can
    be used by the user to create their own FragmentationModels.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (groups names). Mandatory columns:
        'detection_smarts' (SMARTS representations of the group to detect its
        precense in the molecule), 'smarts' (true SMARTS of the group witouht
        additional atom detections), 'contribute' (dictionary as a string with
        the group contribution), 'composed' (y or n if it is or is not a
        composed structure), 'molecular_weight' (molecular weight of the
        subgroups).
    split_detection_smarts : List[str], optional
        List of subgroups that have different SMARTS representations. by
        default []
    problematic_structures : Union[pd.DataFrame, None], optional
        Model's problematic/ambiguous structures. Index: 'smarts' (SMARTS of
        the problematic structure). Mandatory columns: 'contribute' (dictionary
        as a string with the structure contribution), by default None
    hideouts : Union[pd.DataFrame, None], optional
        Hideouts for each group. Index: 'group' (Group of the model that can be
        hiden). Mandatory columns: 'hideout' (other subgroups to find the hiden
        subgroup), by defautl None

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (groups names). Mandatory columns:
        'detection_smarts' (SMARTS representations of the group to detect its
        precense in the molecule), 'smarts' (true SMARTS of the group witouht
        additional atom detections. If a value is missing uses the
        corresponding detection_smarts), 'contribute' (dictionary as a string
        with the group contribution), 'composed' (y or n if it is or is not a
        composed structure), 'molecular_weight' (molecular weight of the
        subgroups).
    split_detection_smarts : List[str]
        List of subgroups that have different SMARTS representations.
    problematic_structures : pd.Dataframe
        Model's problematic/ambiguous structures. Index: 'smarts' (SMARTS of
        the problematic structure). Mandatory columns: 'contribute' (dictionary
        as a string with the structure contribution)
    hideouts : pd.DataFrame
        Hideouts for each group. Index: 'group' (Group of the model that can be
        hiden). Mandatory columns: 'hideout' (other subgroups to find the hiden
        subgroup)
    detection_mols : dict
        Dictionary cotaining all the rdkit Mol object from the detection_smarts
        subgroups column.
    fit_mols : dict
        Dictionary cotaining all the rdkit Mol object from the smarts subgroups
        column.
    contribution_matrix : pd.Dataframe
        Contribution matrix of the model built from the subgroups contribute.
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        split_detection_smarts: List[str] = [],
        problematic_structures: Union[pd.DataFrame, None] = None,
        hideouts: Union[pd.DataFrame, None] = None,
    ) -> None:
        self.subgroups = subgroups
        self.split_detection_smarts = split_detection_smarts

        # =====================================================================
        # Empty problematics template
        # =====================================================================
        if problematic_structures is None:
            self.problematic_structures = pd.DataFrame(
                [], columns=["smarts", "contribute"]
            ).set_index("smarts")
        else:
            self.problematic_structures = problematic_structures

        # =====================================================================
        # Hideouts
        # =====================================================================
        if hideouts is None:
            self.hideouts = pd.DataFrame(
                [], columns=["group", "hideout"]
            ).set_index("group")
        else:
            self.hideouts = hideouts

        # =====================================================================
        # Contribution matrix build
        # =====================================================================
        self.contribution_matrix = self._build_contrib_matrix()

        # =====================================================================
        # Instantiate all de mol object from their smarts
        # =====================================================================
        self.detection_mols = self._instantiate_detection_mol()
        self.fit_mols = self._instantiate_fit_mols()

    def _build_contrib_matrix(self) -> pd.DataFrame:
        """Build the contribution matrix of the model.

        Returns
        -------
        pd.DataFrame
            The contribution matrix of the model built from the contribute
            column of the subgroups DataFrame.

        Raises
        ------
        ValueError
            Bad contribution parsing of a group
        TypeError
            Bad contribution parsing of a group
        """
        index = self.subgroups.index.to_numpy()
        matrix = np.zeros((len(index), len(index)), dtype=int)

        # Build the matrix
        dfm = pd.DataFrame(matrix, index=index, columns=index).rename_axis(
            "group"
        )

        # Fill the matrix
        for group in index:
            str_contribution = self.subgroups.loc[group, "contribute"]

            try:
                contribution = json.loads(str_contribution)
            except json.JSONDecodeError:
                raise ValueError(
                    f"Bad contribute parsing of the group: {group}"
                )
            except TypeError:
                raise TypeError(
                    f"Bad contribution parsing of the group: {group}."
                )

            for k in contribution.keys():
                dfm.loc[group, k] = contribution[k]

        return dfm

    def _instantiate_detection_mol(self) -> dict:
        """Instantiate all the rdkit Mol object from the detection_smarts.

        Returns
        -------
        dict
            Mol objects.
        """
        mols = {}

        for group in self.subgroups.index:
            if group not in self.split_detection_smarts:
                mols[group] = [
                    Chem.MolFromSmarts(
                        self.subgroups.loc[group, "detection_smarts"]
                    )
                ]
            else:
                smarts = self.subgroups.loc[group, "detection_smarts"].split(
                    ","
                )

                mol_smarts = []
                for sms in smarts:
                    mol_smarts.append(Chem.MolFromSmarts(sms))

                mols[group] = mol_smarts
        return mols

    def _instantiate_fit_mols(self) -> dict:
        """Instantiate all the rdkit Mol object from the smarts.

        Returns
        -------
        dict
            Mol object to fit the subgroups in the molecule's atoms.
        """
        mols = {}

        for group in self.subgroups.index:
            smarts = self.subgroups.loc[group, "smarts"]

            if isinstance(smarts, str):
                mols[group] = [Chem.MolFromSmarts(smarts)]
            else:
                mols[group] = self.detection_mols[group]

        return mols
