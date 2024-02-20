"""FragmentationModel module.

All ugropy models (dortmund, joback, unifac, psrk) are instances of the
FragmentationModule class.
"""

import json
from typing import List, Union

import numpy as np

import pandas as pd


class FragmentationModel:
    """FragmentationModel class.

    All ugropy supported models are an instance of this class. This class can
    be used by the user to create their own FragmentationModels.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (groups names). Mandatory columns:
        'smarts' (SMARTS representations of the group), 'contribute'
        (dictionary as a string with the group contribution), 'composed' (y or
        n if it is or is not a composed structure), 'molecular_weight'
        (molecular weight of the subgroups).
    main_groups : Union[pd.DataFrame, None], optional
        List of main groups, by now it is optional and does nothing, by default
        None
    problematic_structures : Union[pd.DataFrame, None], optional
        Model's problematic/ambiguous structures. Index: 'smarts' (SMARTS of
        the problematic structure). Mandatory columns: 'contribute' (dictionary
        as a string with the structure contribution), by default None
    ch2_hideouts : List[str], optional
        _description_, by default []
    ch_hideouts : List[str], optional
        _description_, by default []
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        main_groups: Union[pd.DataFrame, None] = None,
        problematic_structures: Union[pd.DataFrame, None] = None,
        ch2_hideouts: List[str] = [],
        ch_hideouts: List[str] = [],
        hidding_groups: List[str] = [],
    ):
        self.subgroups = subgroups

        # Empty main_groups template
        if main_groups is None:
            self.main_groups = pd.DataFrame(
                [], columns=["no.", "main group name", "subgroups"]
            ).set_index("no.")
        else:
            self.main_groups = main_groups

        # Empty problematic_structures template
        if problematic_structures is None:
            self.problematic_structures = pd.DataFrame(
                [], columns=["smarts", "contribute"]
            ).set_index("smarts")
        else:
            self.problematic_structures = problematic_structures

        # Hideouts
        self.ch2_hideouts = ch2_hideouts
        self.ch_hideouts = ch_hideouts

        # Contribution matrix
        self.contribution_matrix = self._build_contrib_matrix()

        # Shared hideouts
        self.hidding_groups = hidding_groups
        self.shared_hideouts = self._build_shared_hideouts()

    def _build_contrib_matrix(self) -> pd.DataFrame:
        """Build the contribution matrix of the model.

        Returns
        -------
        pd.DataFrame
            The contribution matrix of the model built from the contribute
            column of the subgroups DataFrame.

        Raises
        ------
        json.JSONDecodeError
            Bad contribution parsing of a group
        TypeError
            Bad contribution parsing of a group
        """
        index = self.subgroups.index.to_numpy()
        matrix = np.zeros((len(index), len(index)), dtype=int)

        # build the matrix
        dfm = pd.DataFrame(matrix, index=index, columns=index).rename_axis(
            "group"
        )

        # fill the matrix
        for group in index:
            str_contribution = self.subgroups.loc[group, "contribute"]

            try:
                contribution = json.loads(str_contribution)
            except json.JSONDecodeError:
                raise json.JSONDecodeError(
                    f"Bad contribution parsing of the group: {group}"
                )
            except TypeError:
                raise TypeError(
                    f"Bad contribution parsing of the group: {group}"
                )

            for k in contribution.keys():
                dfm.loc[group, k] = contribution[k]

        return dfm

    def _build_shared_hideouts(self) -> pd.DataFrame:
        shared_hideouts = pd.DataFrame(columns=["shared_group", "hideout"])
        shared_hideouts.set_index("shared_group", inplace=True)

        if len(self.hidding_groups) == 0:
            return shared_hideouts

        subgroups = self.subgroups.index
        contributes = self.subgroups["contribute"]

        for shared in self.hidding_groups:
            for group, contribute in zip(subgroups, contributes):
                dict_contrib = json.loads(contribute)

                if dict_contrib.get(shared, 0) < 0:
                    new_row = pd.DataFrame([(shared, group)], columns=["shared_group", "hideout"])
                    new_row.set_index("shared_group", inplace=True)

                    shared_hideouts = pd.concat([shared_hideouts, new_row])

        return shared_hideouts
