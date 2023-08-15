import ugropy as ug

import pytest


# =============================================================================
# 3- ACH Main group: ACH, AC
# =============================================================================

# UNIFAC
trials_unifac = [
    # phenanthrene
    ("C1=CC=C2C(=C1)C=CC3=CC=CC=C32", {"ACH": 10, "AC": 4}, "smiles"),
    # Anthracene
    ("C1=CC=C2C=C3C=CC=CC3=CC2=C1", {"ACH": 10, "AC": 4}, "smiles"),
    # 1,1-Diphenylethylene
    ("C=C(C1=CC=CC=C1)C2=CC=CC=C2", {"ACH": 10, "AC": 2, "CH2=C": 1}, "smiles"),
    # biphenyl
    ("C1=CC=C(C=C1)C2=CC=CC=C2", {"ACH": 10, "AC": 2}, "smiles"),
    # Naphthalene
    ("C1=CC=C2C=CC=CC2=C1", {"ACH": 8, "AC": 2}, "smiles"),
    ("benzene", {"ACH": 6}, "name"),
    ("styrene", {"AC": 1, "CH2=CH": 1, "ACH": 5}, "name")
]

@pytest.mark.ACH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ach_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result