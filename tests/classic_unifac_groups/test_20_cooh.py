import pytest

import ugropy as ug


# =============================================================================
# 20- COOH Main group: COOH, HCOOH
# =============================================================================

# UNIFAC
trials_unifac = [
    # 2,4-Diaminobutyric acid
    (
        "C(CN)C(C(=O)O)N",
        {"COOH": 1, "CHNH2": 1, "CH2": 1, "CH2NH2": 1},
        "smiles",
    ),
    # acetic acid
    ("CC(=O)O", {"CH3": 1, "COOH": 1}, "smiles"),
    # formic acid
    ("C(=O)O", {"HCOOH": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cooh_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result
