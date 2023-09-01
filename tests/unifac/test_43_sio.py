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
    ("CC(C)O[Si](C)(C)C", {"CH3": 5, "CH-O": 1, "SI": 1}, "smiles"),
    ("CC(C)O[SiH](C)C", {"CH3": 4, "SIH": 1, "CH-O": 1}, "smiles"),
    ("C[SiH2]OC(C)C", {"CH3": 3, "SIH2": 1, "CH-O": 1}, "smiles"),
    ("CCO[SiH2]C", {"CH3": 2, "SIH2": 1, "CH2O": 1}, "smiles"),
    ("CO[SiH2]C", {"CH3": 1, "SIH2": 1, "CH3O": 1}, "smiles"),
    ("C[Si](O[Si](C)([H])[H])([H])[H]", {"CH3": 2, "SIH2O": 1, "SIH2": 1}, "smiles"),
    ("C[Si](C)(O[Si](C)(C)[H])[H]", {"CH3": 4, "SIHO": 1, "SIH": 1}, "smiles"),
    ("Octamethylcyclotetrasiloxane", {"CH3": 8, "SIO": 4}, "name"),
    # Esters + ether
    ("CC(=O)O[SiH3]", {"CH3COO": 1, "SIH3": 1}, "smiles"),
    ("CC(=O)O[SiH2][SiH3]", {"CH3COO": 1, "SIH3": 1, "SIH2": 1}, "smiles"),
    ("CC(=O)O[SiH]([SiH3])[SiH3]", {"CH3COO": 1, "SIH3": 2, "SIH": 1}, "smiles"),
    ("CC(=O)O[Si]([SiH3])([SiH3])[SiH3]", {"CH3COO": 1, "SIH3": 3, "SI": 1}, "smiles"),
    ("CCC(=O)O[Si]([SiH3])([SiH3])[SiH3]", {"CH3": 1, "CH2COO": 1, "SIH3": 3, "SI": 1}, "smiles"),
    ("CC(C)C(=O)O[Si]([SiH3])([SiH3])[SiH3]", {"CH3": 2, "CH": 1, "COO": 1, "SIH3": 3, "SI": 1}, "smiles"),
    ("[SiH3]OC(=O)O[SiH2][SiH3]", {"COO": 1, "SIH3": 2, "SIH2O": 1}, "smiles"),
    ("[SiH3][SiH2]OC(=O)O[SiH2][SiH3]", {"COO": 1, "SIH3": 2, "SIH2": 1, "SIH2O": 1}, "smiles"),
    ("[SiH3][SiH2]OC(=O)O[SiH]([SiH3])[SiH3]", {"COO": 1, "SIH3": 3, "SIH": 1, "SIH2O": 1}, "smiles"),
    ("[SiH3][SiH2]OC(=O)O[Si]([SiH3])([SiH3])[SiH3]", {"COO": 1, "SIH3": 4, "SI": 1, "SIH2O": 1}, "smiles"),
    ("[SiH3]OC(=O)O[SiH]([SiH3])[SiH3]", {"COO": 1, "SIH3": 3, "SIHO": 1}, "smiles"),
    ("[SiH3][SiH]([SiH3])OC(=O)O[SiH]([SiH3])[SiH3]", {"COO": 1, "SIH3": 4, "SIH": 1, "SIHO": 1}, "smiles"),
    ("[SiH3][SiH]([SiH3])OC(=O)O[Si]([SiH3])([SiH3])[SiH3]", {"COO": 1, "SIH3": 5, "SI": 1, "SIHO": 1}, "smiles"),
    ("[SiH3]OC(=O)O[Si]([SiH3])([SiH3])[SiH3]", {"COO": 1, "SIH3": 4, "SIO": 1}, "smiles"),
    ("[SiH3][Si]([SiH3])([SiH3])OC(=O)O[Si]([SiH3])([SiH3])[SiH3]", {"COO": 1, "SIH3": 6, "SI": 1, "SIO": 1}, "smiles"),
    # Ester + ether (mix C and Si)
    # TODO: mixup SI and C
    # I trully hate ether group
    ("[SiH3][SiH]([SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]", {"SIH3": 4, "SIHO": 1, "SIH2O": 1, "SIH": 1}, "smiles"),
    ("[SiH3][SiH]([SiH3])O[SiH2]O[SiH](O[SiH2]O[SiH]([SiH3])[SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]", {"SIH3": 6, "SIH": 1, "SIH2O": 3, "SIHO": 3}, "smiles"),
    ("[SiH3][SiH]([SiH3])O[SiH2]O[SiH2][SiH](O[SiH2]O[SiH]([SiH3])[SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]", {"SIH3": 6, "SIH": 2, "SIH2O": 4, "SIHO": 2}, "smiles"),
    ("[SiH3][SiH]([SiH3])O[SiH2]O[SiH]([SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]", {"SIH3": 5, "SIH": 1, "SIH2O": 2, "SIHO": 2}, "smiles"),
    ("[SiH3]O[SiH2]O[SiH]([SiH3])O[SiH2]O[SiH3]", {}, "smiles"),
    ("[SiH3][SiH](O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])([SiH3])[SiH3]", {"SIH3": 7, "SI": 1, "SIHO": 1, "SIO": 1}, "smiles"),
    ("[SiH3][SiH](O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])(O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3])O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3]", {"SIH3": 13, "SI": 1, "SIHO": 3, "SIO": 3}, "smiles"),
    ("[SiH3][SiH](O[SiH2][Si]([SiH3])(O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3])O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])([SiH3])[SiH3]", {"SIH3": 13, "SI": 2, "SIH2O": 1, "SIHO": 3, "SIO": 2}, "smiles"),
    ("[SiH3][SiH](O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])([SiH3])O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3]", {"SIH3": 10, "SI": 1, "SIHO": 2, "SIO": 2}, "smiles"),
    # Ether concatenation (mix C and Si)
    # TODO: mixup SI and C
]

@pytest.mark.SIO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sio_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result