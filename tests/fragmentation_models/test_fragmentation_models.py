import pandas as pd

import pytest

from ugropy import FragmentationModel, unifac
from ugropy.fragmentation_models import models


def test_no_problem_model():
    unifac_no_problem = FragmentationModel(
        subgroups=models._uni_sg, hideouts=models._uni_hide
    )

    assert unifac_no_problem.subgroups.equals(unifac.subgroups)
    assert unifac_no_problem.contribution_matrix.equals(
        unifac.contribution_matrix
    )

    df = pd.DataFrame([], columns=["smarts", "contribute"]).set_index("smarts")

    assert unifac_no_problem.problematic_structures.equals(df)


def test_making_it_explode():
    subgroups = models._uni_sg.copy()
    subgroups.loc["CH3", "contribute"] = "pitágoras"

    with pytest.raises(ValueError):
        FragmentationModel(subgroups=subgroups)

    subgroups.loc["CH3", "contribute"] = 2

    with pytest.raises(TypeError):
        FragmentationModel(subgroups=subgroups)

    subgroups.loc["CH3", "contribute"] = ""

    with pytest.raises(ValueError):
        FragmentationModel(subgroups=subgroups)
