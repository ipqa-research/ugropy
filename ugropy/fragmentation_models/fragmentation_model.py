"""FragmentationModel module.

All ugropy models (dortmund, joback, unifac, psrk) are instances of the
FragmentationModule class.
"""

import json
from typing import Union

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
    hideouts : Union[pd.DataFrame, None], optional
        Hideouts for each group. Mandatory columns: 'group' (Group of the model
        that can be hiden), 'hideout' (other subgroups to find the hiden
        subgroup), by default []

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (groups names). Mandatory columns:
        'smarts' (SMARTS representations of the group), 'contribute'
        (dictionary as a string with the group contribution), 'composed' (y or
        n if it is or is not a composed structure), 'molecular_weight'
        (molecular weight of the subgroups).
    main_groups : Union[pd.DataFrame, None]
        List of main groups, by now it is optional and does nothing
    problematic_structures : Union[pd.DataFrame, None]
        Model's problematic/ambiguous structures. Index: 'smarts' (SMARTS of
        the problematic structure). Mandatory columns: 'contribute' (dictionary
        as a string with the structure contribution)
    hideouts : Union[pd.DataFrame, None]
        Hideouts for each group. Mandatory columns: 'group' (Group of the model
        that can be hiden), 'hideout' (other subgroups to find the hiden
        subgroup)
    contribution_matrix : pd.DataFrame
        Model's contribution matrix built from subgroups contribute.
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        main_groups: Union[pd.DataFrame, None] = None,
        problematic_structures: Union[pd.DataFrame, None] = None,
        hideouts: Union[pd.DataFrame, None] = None,
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
        if hideouts is None:
            self.hideouts = pd.DataFrame(
                [], columns=["group", "hideout"]
            ).set_index("group")
        else:
            self.hideouts = hideouts

        # Contribution matrix
        self.contribution_matrix = self._build_contrib_matrix()

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
