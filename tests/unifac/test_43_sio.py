import ugropy as ug

import pytest


# =============================================================================
# 43- SIO Main group: SIH2O, SIHO, SIO
# =============================================================================

# UNIFAC
trials_unifac = [
    ("C[SiH2]O[SiH](C)C",  {"CH3": 3, "SIH2O": 1, "SIH": 1}, "smiles"),
    ("C[SiH2]O[Si](C)(C)C", {"CH3": 4, "SIH2O": 1, "SI": 1}, "smiles"),
    ("C[SiH](C)O[Si](C)(C)C", {"CH3": 5, "SIHO": 1, "SI": 1}, "smiles"),
    ("CC(C)(C)O[Si](C)(C)C", {"CH3": 6, "SIO": 1, "C": 1}, "smiles"),
    ("CC(C)O[Si](C)(C)C", {"CH3": 5, "SIO": 1, "CH": 1}, "smiles"),
    ("CC(C)O[SiH](C)C", {"CH3": 4, "SIHO": 1, "CH": 1}, "smiles"),
    ("C[SiH2]OC(C)C", {"CH3": 3, "SIH2O": 1, "CH": 1}, "smiles"),
    ("CCO[SiH2]C", {"CH3": 2, "SIH2O": 1, "CH2": 1}, "smiles"),
    ("CO[SiH2]C", {"CH3O": 1, "SIH2": 1, "CH3": 1}, "smiles"),
    ("C[Si](O[Si](C)([H])[H])([H])[H]", {"CH3": 2, "SIH2O": 1, "SIH2": 1}, "smiles"),
    ("C[Si](C)(O[Si](C)(C)[H])[H]", {"CH3": 4, "SIHO": 1, "SIH": 1}, "smiles"),
    ("Octamethylcyclotetrasiloxane", {"CH3": 8, "SIO": 4}, "name"),
]

@pytest.mark.SIO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sio_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result