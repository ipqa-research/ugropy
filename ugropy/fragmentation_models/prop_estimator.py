"""PropertiesEstimator module."""

from typing import List, Union

import pandas as pd

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


class PropertiesEstimator(FragmentationModel):
    """Fragmentation model dedicated to properties estimation models.

    joback is a instance of this class.

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
    properties_contributions : pd.DataFrame, optional
        Group's properties contributions, by default None

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
    properties_contributions : pd.DataFrame
        Group's properties contributions.
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        split_detection_smarts: List[str] = [],
        problematic_structures: Union[pd.DataFrame, None] = None,
        hideouts: Union[pd.DataFrame, None] = None,
        properties_contributions: Union[pd.DataFrame, None] = None,
    ):
        super().__init__(
            subgroups, split_detection_smarts, problematic_structures, hideouts
        )

        # =====================================================================
        # Properties contributions
        # =====================================================================
        self.properties_contributions = properties_contributions
