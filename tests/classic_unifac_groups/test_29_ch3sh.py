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


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3sh_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result

    if identifier != "CCCC1=CC=C(C=C1)C(C)S":
        assert ug.get_groups(ug.psrk, identifier, identifier_type) == result
    else:
        ug.get_psrk_groups(identifier, identifier_type) == {
            "CH3": 2,
            "CH2": 1,
            "ACH": 4,
            "ACCH2": 1,
            "CHSH": 1,
            "AC": 1,
        }
