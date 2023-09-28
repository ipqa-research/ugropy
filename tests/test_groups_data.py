import re
from pathlib import Path

import pandas as pd

from ugropy import constants


def test_unifac_data():
    here = Path(__file__).parent.resolve()

    with open(f"{here}/original_unifac.csv", mode="r") as f:
        df = pd.read_csv(f, sep="|", index_col="Subgroup Name", comment="?")

    for group in df.index:
        assert df.loc[group, "R"] == constants.unifac_subgroups.loc[group, "R"]
        assert df.loc[group, "Q"] == constants.unifac_subgroups.loc[group, "Q"]
        assert (
            df.loc[group, "No."]
            == constants.unifac_subgroups.loc[group, "subgroup_number"]
        )

        pattern = r"\[(\d+)\]"
        main_group = df.loc[group, "Maingroup"]
        main_group_num = int(re.search(pattern, main_group).group(1))
        assert (
            main_group_num
            == constants.unifac_subgroups.loc[group, "main_group"]
        )


def test_psrk_data():
    here = Path(__file__).parent.resolve()

    with open(f"{here}/original_psrk.csv", mode="r") as f:
        df = pd.read_csv(f, sep="|", index_col="Subgroup name", comment="?")

    for group in df.index:
        assert df.loc[group, "R"] == constants.psrk_subgroups.loc[group, "R"]
        assert df.loc[group, "Q"] == constants.psrk_subgroups.loc[group, "Q"]
        assert (
            df.loc[group, "No."]
            == constants.psrk_subgroups.loc[group, "subgroup_number"]
        )

        pattern = r"\[(\d+)\]"
        main_group = df.loc[group, "Main group"]
        main_group_num = int(re.search(pattern, main_group).group(1))
        assert (
            main_group_num == constants.psrk_subgroups.loc[group, "main_group"]
        )


def test_dormund_data():
    here = Path(__file__).parent.resolve()

    with open(f"{here}/original_dortmund.csv", mode="r") as f:
        df = pd.read_csv(f, sep="|", index_col="Subgroup Name", comment="?")

    for group in df.index:
        assert (
            df.loc[group, "R"] == constants.dortmund_subgroups.loc[group, "R"]
        )
        assert (
            df.loc[group, "Q"] == constants.dortmund_subgroups.loc[group, "Q"]
        )
        assert (
            df.loc[group, "No."]
            == constants.dortmund_subgroups.loc[group, "subgroup_number"]
        )
        assert (
            df.loc[group, "Main Group No."]
            == constants.dortmund_subgroups.loc[group, "main_group"]
        )
