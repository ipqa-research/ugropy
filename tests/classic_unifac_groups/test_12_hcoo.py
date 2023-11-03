import pytest

import ugropy as ug


# =============================================================================
# 12- HCOO Main group: HCOO
# =============================================================================

# UNIFAC
trials_unifac = [
    # phenyl formate
    ("C1=CC=C(C=C1)OC=O", {"ACH": 5, "AC": 1, "HCOO": 1}, "smiles"),
    # methyl formate
    ("COC=O", {"HCOO": 1, "CH3": 1}, "smiles"),
    # ethyl formate
    ("CCOC=O", {"HCOO": 1, "CH2": 1, "CH3": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cho_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result
