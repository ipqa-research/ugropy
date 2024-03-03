from pathlib import Path

import pandas as pd

import pytest

from ugropy import FragmentationModel, get_groups, unifac
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


def test_draw():
    here = Path(__file__).parent.resolve()

    with open(f"{here}/../../logo.svg") as f:
        expected = f.read()

    mol = get_groups(
        unifac, "CCCC1=C(COC(C)(C)COC(=O)OCC)C=C(CC2=CC=CC=C2)C=C1", "smiles"
    )

    svg = mol.draw("ugropy", 800, 450, 50, 14)

    assert expected == svg
