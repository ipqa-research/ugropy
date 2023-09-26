import pytest

import ugropy as ug


# =============================================================================
# 29- CH3SH Main group: CH3SH, CH2SH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("CCCC1=CC=C(C=C1)C(C)S", {}, "smiles"),
    (
        "CCCC1=CC=C(CS)C=C1",
        {"ACCH2": 1, "ACH": 4, "AC": 1, "CH2SH": 1, "CH2": 1, "CH3": 1},
        "smiles",
    ),
    (
        "CC1=CC=C(CS)C=C1",
        {"ACCH3": 1, "ACH": 4, "AC": 1, "CH2SH": 1},
        "smiles",
    ),
    ("SCC1=CC=NC=C1", {"C5H4N": 1, "CH2SH": 1}, "smiles"),
    # methanethiol
    ("CS", {"CH3SH": 1}, "smiles"),
    # ethanethiol
    ("CCS", {"CH2SH": 1, "CH3": 1}, "smiles"),
    # impossible
    ("CCCC1=CC=C(C=C1)C(C)S", {}, "smiles"),
]


@pytest.mark.CH3SH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3sh_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result

    if identifier != "CCCC1=CC=C(C=C1)C(C)S":
        assert groups.psrk_groups == result
    else:
        assert groups.psrk_groups == {
            "CH3": 2,
            "CH2": 1,
            "ACH": 4,
            "ACCH2": 1,
            "CHSH": 1,
            "AC": 1,
        }
