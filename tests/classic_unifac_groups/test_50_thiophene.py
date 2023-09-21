import pytest

import ugropy as ug


# =============================================================================
# 50- thiophene Main group: C4H4S, C4H3S, C4H2S
# =============================================================================

# UNIFAC
trials_unifac = [
    ("thiophene", {"C4H4S": 1}, "name"),
    ("2-methylthiophene", {"C4H3S": 1, "CH3": 1}, "name"),
    ("2,3-dimethylthiophene", {"C4H2S": 1, "CH3": 2}, "name"),
    ("OC1=CSC=C1", {"C4H3S": 1, "OH": 1}, "smiles"),
    ("OC1=CC=CS1", {"C4H3S": 1, "OH": 1}, "smiles"),
    ("OC1=CSC=C1O", {"C4H2S": 1, "OH": 2}, "smiles"),
    ("OC1=CC(O)=CS1", {"C4H2S": 1, "OH": 2}, "smiles"),
    ("OC1=CC=C(O)S1", {"C4H2S": 1, "OH": 2}, "smiles"),
]


@pytest.mark.thiophene
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_thiophene_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
