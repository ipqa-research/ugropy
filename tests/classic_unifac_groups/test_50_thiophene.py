import pytest

import ugropy as ug


# =============================================================================
# 50- thiophene Main group: C4H4S, C4H3S, C4H2S
# =============================================================================

# UNIFAC
trials_unifac = [
    ("CC1=C(SC=C1)C", {"C4H2S": 1, "CH3": 2}, "smiles"),
    # thiophene
    ("C1=CSC=C1", {"C4H4S": 1}, "smiles"),
    # 2-methylthiophene
    ("CC1=CC=CS1", {"C4H3S": 1, "CH3": 1}, "smiles"),
    # 2,3-dimethylthiophene
    ("OC1=CSC=C1", {"C4H3S": 1, "OH": 1}, "smiles"),
    ("OC1=CC=CS1", {"C4H3S": 1, "OH": 1}, "smiles"),
    ("OC1=CSC=C1O", {"C4H2S": 1, "OH": 2}, "smiles"),
    ("OC1=CC(O)=CS1", {"C4H2S": 1, "OH": 2}, "smiles"),
    ("OC1=CC=C(O)S1", {"C4H2S": 1, "OH": 2}, "smiles"),
]


@pytest.mark.PSRK
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_thiophene_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result
