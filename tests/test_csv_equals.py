from ugropy import constants


def test_unifac_psrk_equal():
    """check if the same subgroup appears in 2 models, both xxx_subgroups.csv
    rows must be the same.
    """
    for group in constants.unifac_subgroups.index:
        if group in constants.psrk_subgroups.index:
            assert (
                constants.unifac_subgroups.loc[group, "smarts"]
                == constants.psrk_subgroups.loc[group, "smarts"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "contribute"]
                == constants.psrk_subgroups.loc[group, "contribute"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "composed"]
                == constants.psrk_subgroups.loc[group, "composed"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "molecular_weight"]
                == constants.psrk_subgroups.loc[group, "molecular_weight"]
            )


def test_unifac_dortmund_equal():
    """check if the same subgroup appears in 2 models, both xxx_subgroups.csv
    rows must be the same.
    """
    for group in constants.unifac_subgroups.index:
        if group in constants.dortmund_subgroups.index:
            assert (
                constants.unifac_subgroups.loc[group, "smarts"]
                == constants.dortmund_subgroups.loc[group, "smarts"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "contribute"]
                == constants.dortmund_subgroups.loc[group, "contribute"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "composed"]
                == constants.dortmund_subgroups.loc[group, "composed"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "molecular_weight"]
                == constants.dortmund_subgroups.loc[group, "molecular_weight"]
            )


def test_psrk_dortmund_equal():
    """check if the same subgroup appears in 2 models, both xxx_subgroups.csv
    rows must be the same.
    """
    for group in constants.psrk_subgroups.index:
        if group in constants.dortmund_subgroups.index:
            assert (
                constants.psrk_subgroups.loc[group, "smarts"]
                == constants.dortmund_subgroups.loc[group, "smarts"]
            )

            assert (
                constants.psrk_subgroups.loc[group, "contribute"]
                == constants.dortmund_subgroups.loc[group, "contribute"]
            )

            assert (
                constants.psrk_subgroups.loc[group, "composed"]
                == constants.dortmund_subgroups.loc[group, "composed"]
            )

            assert (
                constants.psrk_subgroups.loc[group, "molecular_weight"]
                == constants.dortmund_subgroups.loc[group, "molecular_weight"]
            )
