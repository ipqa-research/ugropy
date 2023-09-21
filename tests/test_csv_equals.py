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
                constants.unifac_subgroups.loc[group, "subgroup_number"]
                == constants.psrk_subgroups.loc[group, "subgroup_number"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "main_group"]
                == constants.psrk_subgroups.loc[group, "main_group"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "R"]
                == constants.psrk_subgroups.loc[group, "R"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "Q"]
                == constants.psrk_subgroups.loc[group, "Q"]
            )

            assert (
                constants.unifac_subgroups.loc[group, "molecular_weight"]
                == constants.psrk_subgroups.loc[group, "molecular_weight"]
            )
