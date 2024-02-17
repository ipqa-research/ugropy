import pytest

import ugropy as ug


# =============================================================================
# 29|CH3SH|[59]CH3SH [60]CH2SH
# =============================================================================

# Dortmund
trials = [
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
    ("SCC1=CC=NC=C1", {"ACH": 2, "AC2H2N": 1, "CH2SH": 1, "AC": 1}, "smiles"),
    # methanethiol
    ("CS", {"CH3SH": 1}, "smiles"),
    # ethanethiol
    ("CCS", {"CH2SH": 1, "CH3": 1}, "smiles"),
    # impossible
    ("CCCC1=CC=C(C=C1)C(C)S", {}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ch3sh_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result
