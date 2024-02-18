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

    All ugropy FragmentationModels are an instance of this class. This class
    can be used by the user to create their own FragmentationModels.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups.
    main_groups : Union[pd.DataFrame, None], optional
        _description_, by default None
    problematic_structures : Union[pd.DataFrame, None], optional
        _description_, by default None
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

    def _build_contrib_matrix(self) -> pd.DataFrame:
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
                    f"Band contribution parsing of the group: {group}"
                )
            except TypeError:
                raise json.JSONDecodeError(
                    f"Band contribution parsing of the group: {group}"
                )

            for k in contribution.keys():
                dfm.loc[group, k] = contribution[k]

        return dfm
