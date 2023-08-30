import ugropy as ug

import pytest


# =============================================================================
# 42- SIH2 Main group: SIH3, SIH2, SIH, SI
# =============================================================================

# UNIFAC
trials_unifac = [
    ("methylsilane", {"CH3": 1, "SIH3": 1}, "name"),
    ("CC[Si](CC)([H])[H]", {"CH3": 2, "CH2": 2, "SIH2": 1}, "smiles"),
    ("C[Si](O[Si](C)(C)C)(O[Si](C)(C)C)[H]", {"CH3": 7, "SIO": 1, "SIHO": 1, "SI": 1}, "smiles"),
    ("Hexamethyldisiloxane", {"CH3": 6, "SIO": 1, "SI": 1}, "name"),
]

@pytest.mark.SIH2
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sih2_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result