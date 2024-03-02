import pandas as pd

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


class GibbsModel(FragmentationModel):
    def __init__(
        self,
        subgroups: pd.DataFrame,
        main_groups: pd.DataFrame | None = None,
        problematic_structures: pd.DataFrame | None = None,
        hideouts: pd.DataFrame | None = None,
    ):
        super().__init__(
            subgroups, main_groups, problematic_structures, hideouts
        )

        # =====================================================================
        # Empty main_groups DataFrame template
        # =====================================================================
        if main_groups is None:
            self.main_groups = pd.DataFrame(
                [], columns=["no.", "main group name", "subgroups"]
            ).set_index("no.")
        else:
            self.main_groups = main_groups
