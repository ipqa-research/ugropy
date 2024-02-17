import json

from ugropy import constants


# test if i'm not missing any potential hideout
# =============================================================================
# UNIFAC
# =============================================================================
def test_unifac_ch2_hideouts():
    for group in constants.unifac_subgroups.index:
        # exceptions (full molecules gropus and well ch2...)
        if group in ["CH2", "THF", "DOH", "NMP", "MORPH"]:
            continue

        string_contribution = constants.unifac_subgroups.loc[
            group, "contribute"
        ]
        contribution = json.loads(string_contribution)

        # Check if CH2 is being substracted in the contribution of the group
        if "CH2" in contribution.keys():
            assert group in constants.unifac_ch2_hide
            continue
        else:
            # Check if CH2 is being substracted in the contribution of the
            # groups in the contribution of "group"
            for grp in list(contribution.keys())[1:]:
                str_contribution = constants.unifac_subgroups.loc[
                    grp, "contribute"
                ]
                contrib = json.loads(str_contribution)

                if "CH2" in contrib.keys():
                    assert group in constants.unifac_ch2_hide
                    break


def test_unifac_ch_hideouts():
    for group in constants.unifac_subgroups.index:
        # exceptions (full molecules gropus and well ch...)
        if group in ["CH", "FURFURAL"]:
            continue

        string_contribution = constants.unifac_subgroups.loc[
            group, "contribute"
        ]
        contribution = json.loads(string_contribution)

        # Check if CH is being substracted in the contribution of the group
        if "CH" in contribution.keys():
            assert group in constants.unifac_ch_hide
            continue
        else:
            # Check if CH is being substracted in the contribution of the
            # groups in the contribution of "group"
            for grp in list(contribution.keys())[1:]:
                str_contribution = constants.unifac_subgroups.loc[
                    grp, "contribute"
                ]
                contrib = json.loads(str_contribution)

                if "CH" in contrib.keys():
                    assert group in constants.unifac_ch_hide
                    break


# =============================================================================
# PSRK
# =============================================================================
def test_psrk_ch2_hideouts():
    for group in constants.psrk_subgroups.index:
        # exceptions (full molecules gropus and well ch2...)
        if group in ["CH2", "THF", "DOH", "NMP", "MORPH", "H2COCH2"]:
            continue

        string_contribution = constants.psrk_subgroups.loc[group, "contribute"]
        contribution = json.loads(string_contribution)

        # Check if CH2 is being substracted in the contribution of the group
        if "CH2" in contribution.keys():
            assert group in constants.psrk_ch2_hide
            continue
        else:
            # Check if CH2 is being substracted in the contribution of the
            # groups in the contribution of "group"
            for grp in list(contribution.keys())[1:]:
                str_contribution = constants.psrk_subgroups.loc[
                    grp, "contribute"
                ]
                contrib = json.loads(str_contribution)

                if "CH2" in contrib.keys():
                    assert group in constants.psrk_ch2_hide
                    break


def test_psrk_ch_hideouts():
    for group in constants.psrk_subgroups.index:
        # exceptions (full molecules gropus and well ch...)
        if group in ["CH", "FURFURAL"]:
            continue

        string_contribution = constants.psrk_subgroups.loc[group, "contribute"]
        contribution = json.loads(string_contribution)

        # Check if CH is being substracted in the contribution of the group
        if "CH" in contribution.keys():
            assert group in constants.psrk_ch_hide
            continue
        else:
            # Check if CH is being substracted in the contribution of the
            # groups in the contribution of "group"
            for grp in list(contribution.keys())[1:]:
                str_contribution = constants.psrk_subgroups.loc[
                    grp, "contribute"
                ]
                contrib = json.loads(str_contribution)

                if "CH" in contrib.keys():
                    assert group in constants.psrk_ch_hide
                    break
