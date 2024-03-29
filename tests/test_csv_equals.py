from ugropy import psrk, unifac


def test_unifac_psrk_equal():
    """check if the same subgroup appears in 2 models, both xxx_subgroups.csv
    rows must be the same.
    """
    for group in unifac.subgroups.index:
        if group in psrk.subgroups.index:
            assert (
                unifac.subgroups.loc[group, "detection_smarts"]
                == psrk.subgroups.loc[group, "detection_smarts"]
            )

            assert (
                unifac.subgroups.loc[group, "contribute"]
                == psrk.subgroups.loc[group, "contribute"]
            )

            assert (
                unifac.subgroups.loc[group, "composed"]
                == psrk.subgroups.loc[group, "composed"]
            )

            assert (
                unifac.subgroups.loc[group, "molecular_weight"]
                == psrk.subgroups.loc[group, "molecular_weight"]
            )

            # assert (
            #     unifac.subgroups["true_smarts"]
            #     .dropna()
            #     .equals(psrk.subgroups["true_smarts"].dropna())
            # )
