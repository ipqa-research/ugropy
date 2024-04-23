import json

from ugropy import psrk, unifac


# test if i'm not missing any potential hideout
# =============================================================================
# UNIFAC
# =============================================================================
def test_unifac_ch2_hideouts():
    for group in unifac.subgroups.index:
        # exceptions (full molecules gropus and well ch2...)
        if group in ["CH2CL2", "CH2", "THF", "DOH", "NMP", "MORPH"]:
            continue

        string_contribution = unifac.subgroups.loc[group, "contribute"]
        contribution = json.loads(string_contribution)

        # Check if CH2 is being substracted in the contribution of the group
        if "CH2" in contribution.keys():
            assert group in unifac.hideouts.loc["CH2"].values
            continue
        else:
            # Check if CH2 is being substracted in the contribution of the
            # groups in the contribution of "group"
            for grp in list(contribution.keys())[1:]:
                str_contribution = unifac.subgroups.loc[grp, "contribute"]
                contrib = json.loads(str_contribution)

                if "CH2" in contrib.keys():
                    assert group in unifac.hideouts.loc["CH2"].values
                    break


def test_unifac_ch_hideouts():
    for group in unifac.subgroups.index:
        # exceptions (full molecules gropus and well ch...)
        if group in ["HCCL2F", "CHCL3", "CH", "FURFURAL"]:
            continue

        string_contribution = unifac.subgroups.loc[group, "contribute"]
        contribution = json.loads(string_contribution)

        # Check if CH is being substracted in the contribution of the group
        if "CH" in contribution.keys():
            assert group in unifac.hideouts.loc["CH"].values
            continue
        else:
            # Check if CH is being substracted in the contribution of the
            # groups in the contribution of "group"
            for grp in list(contribution.keys())[1:]:
                str_contribution = unifac.subgroups.loc[grp, "contribute"]
                contrib = json.loads(str_contribution)

                if "CH" in contrib.keys():
                    assert group in unifac.hideouts.loc["CH"].values
                    break


# =============================================================================
# PSRK
# =============================================================================
def test_psrk_ch2_hideouts():
    for group in psrk.subgroups.index:
        # exceptions (full molecules gropus and well ch2...)
        if group in ["CH2CL2", "CH2", "THF", "DOH", "NMP", "MORPH", "H2COCH2"]:
            continue

        string_contribution = psrk.subgroups.loc[group, "contribute"]
        contribution = json.loads(string_contribution)

        # Check if CH2 is being substracted in the contribution of the group
        if "CH2" in contribution.keys():
            assert group in psrk.hideouts.loc["CH2"].values
            continue
        else:
            # Check if CH2 is being substracted in the contribution of the
            # groups in the contribution of "group"
            for grp in list(contribution.keys())[1:]:
                str_contribution = psrk.subgroups.loc[grp, "contribute"]
                contrib = json.loads(str_contribution)

                if "CH2" in contrib.keys():
                    assert group in psrk.hideouts.loc["CH2"].values
                    break


def test_psrk_ch_hideouts():
    for group in psrk.subgroups.index:
        # exceptions (full molecules gropus and well ch...)
        if group in ["HCCL2F", "CHCL3", "CH", "FURFURAL"]:
            continue

        string_contribution = psrk.subgroups.loc[group, "contribute"]
        contribution = json.loads(string_contribution)

        # Check if CH is being substracted in the contribution of the group
        if "CH" in contribution.keys():
            assert group in psrk.hideouts.loc["CH"].values
            continue
        else:
            # Check if CH is being substracted in the contribution of the
            # groups in the contribution of "group"
            for grp in list(contribution.keys())[1:]:
                str_contribution = psrk.subgroups.loc[grp, "contribute"]
                contrib = json.loads(str_contribution)

                if "CH" in contrib.keys():
                    assert group in psrk.hideouts.loc["CH"].values
                    break
