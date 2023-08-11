import ugropy as ug

import pytest


# =============================================================================
# 43- SIO Main group: SIH2O, SIHO, SIO
# =============================================================================

# UNIFAC
trials_unifac = [
    ("1,3-Dimethyldisiloxane", {"CH3": 2, "SIH2O": 1, "SIH2": 1}, "name"),
    ("1,1,3,3-Tetramethyldisiloxane", {"CH3": 4, "SIHO": 1, "SIH": 1}, "name"),
    ("Octamethylcyclotetrasiloxane", {"CH3": 8, "SIO": 4}, "name"),
    #("C[Si](O[Si](C)([H])[H])([H])[H]", {"CH3": 2, "SIH2O": 1, "SIH2": 1}, "smiles"),
    #("C[Si](C)(O[Si](C)(C)[H])[H]", {"CH3": 4, "SIHO": 1, "SIH": 1}, "smiles"),
    #("Octamethylcyclotetrasiloxane", {"CH3": 8, "SIO": 4}, "name"),
]

@pytest.mark.SIO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sio_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result