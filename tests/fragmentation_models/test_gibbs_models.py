from pathlib import Path

import pandas as pd

from ugropy import unifac
from ugropy.fragmentation_models import GibbsModel


def test_no_info():
    model = GibbsModel(subgroups=unifac.subgroups.copy())

    main_groups = pd.DataFrame(
        [], columns=["no.", "main group name", "subgroups"]
    ).set_index("no.")

    assert main_groups.equals(model.main_groups)

    subgroups_info = pd.DataFrame(
        [],
        columns=["group", "subgroup_number", "main_group", "R", "Q"],
    ).set_index("group")

    assert subgroups_info.equals(model.subgroups_info)
