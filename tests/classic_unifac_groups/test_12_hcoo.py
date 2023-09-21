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
    ("ethyl formate", {"HCOO": 1, "CH2": 1, "CH3": 1}, "name"),
]


@pytest.mark.HCOO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cho_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
