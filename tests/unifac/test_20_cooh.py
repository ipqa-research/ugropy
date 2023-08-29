import ugropy as ug

import pytest


# =============================================================================
# 20- COOH Main group: COOH, HCOOH
# =============================================================================

# UNIFAC
trials_unifac = [
    # 2,4-Diaminobutyric acid
    ("C(CN)C(C(=O)O)N", {"COOH": 1, "CHNH2": 1, "CH2": 1, "CH2NH2": 1}, "smiles"),
    ("acetic acid", {"CH3": 1, "COOH": 1}, "name"),
    ("formic acid", {"HCOOH": 1}, "name"),
]

@pytest.mark.COOH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cooh_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result