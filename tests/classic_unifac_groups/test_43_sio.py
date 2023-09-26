import pytest

import ugropy as ug


# =============================================================================
# 43- SIO Main group: SIH2O, SIHO, SIO
# =============================================================================

# UNIFAC
trials_unifac = [
    ("C[SiH2]O[SiH](C)C", {"CH3": 3, "SIH2O": 1, "SIH": 1}, "smiles"),
    ("C[SiH2]O[Si](C)(C)C", {"CH3": 4, "SIH2O": 1, "SI": 1}, "smiles"),
    ("C[SiH](C)O[Si](C)(C)C", {"CH3": 5, "SIHO": 1, "SI": 1}, "smiles"),
    ("CC(C)(C)O[Si](C)(C)C", {"CH3": 6, "SIO": 1, "C": 1}, "smiles"),
    ("CC(C)O[Si](C)(C)C", {"CH3": 5, "CH-O": 1, "SI": 1}, "smiles"),
    ("CC(C)O[SiH](C)C", {"CH3": 4, "SIH": 1, "CH-O": 1}, "smiles"),
    ("C[SiH2]OC(C)C", {"CH3": 3, "SIH2": 1, "CH-O": 1}, "smiles"),
    ("CCO[SiH2]C", {"CH3": 2, "SIH2": 1, "CH2O": 1}, "smiles"),
    ("CO[SiH2]C", {"CH3": 1, "SIH2": 1, "CH3O": 1}, "smiles"),
    (
        "C[Si](O[Si](C)([H])[H])([H])[H]",
        {"CH3": 2, "SIH2O": 1, "SIH2": 1},
        "smiles",
    ),
    ("C[Si](C)(O[Si](C)(C)[H])[H]", {"CH3": 4, "SIHO": 1, "SIH": 1}, "smiles"),
    # Octamethylcyclotetrasiloxane
    (
        "C[Si]1(O[Si](O[Si](O[Si](O1)(C)C)(C)C)(C)C)C",
        {"CH3": 8, "SIO": 4},
        "smiles",
    ),
    # Esters + ether
    ("CC(=O)O[SiH3]", {"CH3COO": 1, "SIH3": 1}, "smiles"),
    ("CC(=O)O[SiH2][SiH3]", {"CH3COO": 1, "SIH3": 1, "SIH2": 1}, "smiles"),
    (
        "CC(=O)O[SiH]([SiH3])[SiH3]",
        {"CH3COO": 1, "SIH3": 2, "SIH": 1},
        "smiles",
    ),
    (
        "CC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        {"CH3COO": 1, "SIH3": 3, "SI": 1},
        "smiles",
    ),
    (
        "CCC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        {"CH3": 1, "CH2COO": 1, "SIH3": 3, "SI": 1},
        "smiles",
    ),
    (
        "CC(C)C(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        {"CH3": 2, "CH": 1, "COO": 1, "SIH3": 3, "SI": 1},
        "smiles",
    ),
    ("[SiH3]OC(=O)O[SiH2][SiH3]", {"COO": 1, "SIH3": 2, "SIH2O": 1}, "smiles"),
    (
        "[SiH3][SiH2]OC(=O)O[SiH2][SiH3]",
        {"COO": 1, "SIH3": 2, "SIH2": 1, "SIH2O": 1},
        "smiles",
    ),
    (
        "[SiH3][SiH2]OC(=O)O[SiH]([SiH3])[SiH3]",
        {"COO": 1, "SIH3": 3, "SIH": 1, "SIH2O": 1},
        "smiles",
    ),
    (
        "[SiH3][SiH2]OC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        {"COO": 1, "SIH3": 4, "SI": 1, "SIH2O": 1},
        "smiles",
    ),
    (
        "[SiH3]OC(=O)O[SiH]([SiH3])[SiH3]",
        {"COO": 1, "SIH3": 3, "SIHO": 1},
        "smiles",
    ),
    (
        "[SiH3][SiH]([SiH3])OC(=O)O[SiH]([SiH3])[SiH3]",
        {"COO": 1, "SIH3": 4, "SIH": 1, "SIHO": 1},
        "smiles",
    ),
    (
        "[SiH3][SiH]([SiH3])OC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        {"COO": 1, "SIH3": 5, "SI": 1, "SIHO": 1},
        "smiles",
    ),
    (
        "[SiH3]OC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        {"COO": 1, "SIH3": 4, "SIO": 1},
        "smiles",
    ),
    (
        "[SiH3][Si]([SiH3])([SiH3])OC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        {"COO": 1, "SIH3": 6, "SI": 1, "SIO": 1},
        "smiles",
    ),
    ("[SiH3]OC=O", {"HCOO": 1, "SIH3": 1}, "smiles"),
    ("[SiH3][SiH2]OC=O", {"HCOO": 1, "SIH3": 1, "SIH2": 1}, "smiles"),
    ("[SiH3][SiH]([SiH3])OC=O", {"HCOO": 1, "SIH3": 2, "SIH": 1}, "smiles"),
    (
        "[SiH3][Si]([SiH3])([SiH3])OC=O",
        {"HCOO": 1, "SIH3": 3, "SI": 1},
        "smiles",
    ),
    # Ester + ether (mix C and Si)
    ("COC(=O)O[SiH3]", {"CH3O": 1, "COO": 1, "SIH3": 1}, "smiles"),
    ("COC(=O)O[SiH2]C", {"CH3O": 1, "COO": 1, "SIH2": 1, "CH3": 1}, "smiles"),
    ("COC(=O)O[SiH](C)C", {"CH3O": 1, "COO": 1, "SIH": 1, "CH3": 2}, "smiles"),
    (
        "COC(=O)O[Si](C)(C)C",
        {"CH3O": 1, "COO": 1, "SI": 1, "CH3": 3},
        "smiles",
    ),
    ("CCOC(=O)O[SiH3]", {"CH2O": 1, "COO": 1, "SIH3": 1, "CH3": 1}, "smiles"),
    ("CCOC(=O)O[SiH2]C", {"CH2O": 1, "COO": 1, "SIH2": 1, "CH3": 2}, "smiles"),
    (
        "CCOC(=O)O[SiH](C)C",
        {"CH2O": 1, "COO": 1, "SIH": 1, "CH3": 3},
        "smiles",
    ),
    (
        "CCOC(=O)O[Si](C)(C)C",
        {"CH2O": 1, "COO": 1, "SI": 1, "CH3": 4},
        "smiles",
    ),
    (
        "CC(C)OC(=O)O[SiH3]",
        {"CH-O": 1, "COO": 1, "SIH3": 1, "CH3": 2},
        "smiles",
    ),
    (
        "C[SiH2]OC(=O)OC(C)C",
        {"CH-O": 1, "COO": 1, "SIH2": 1, "CH3": 3},
        "smiles",
    ),
    (
        "CC(C)OC(=O)O[SiH](C)C",
        {"CH-O": 1, "COO": 1, "SIH": 1, "CH3": 4},
        "smiles",
    ),
    (
        "CC(C)OC(=O)O[Si](C)(C)C",
        {"CH-O": 1, "COO": 1, "SI": 1, "CH3": 5},
        "smiles",
    ),
    (
        "C[SiH2]OC(=O)OC(C)(C)C",
        {"SIH2O": 1, "COO": 1, "C": 1, "CH3": 4},
        "smiles",
    ),
    (
        "C[SiH](C)OC(=O)OC(C)(C)C",
        {"SIHO": 1, "COO": 1, "C": 1, "CH3": 5},
        "smiles",
    ),
    (
        "CC(C)(C)OC(=O)O[Si](C)(C)C",
        {"SIO": 1, "COO": 1, "C": 1, "CH3": 6},
        "smiles",
    ),
    ("CC(C)(C)OC(=O)O[SiH3]", {}, "smiles"),
    # I trully hate ether group
    (
        "[SiH3][SiH]([SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]",
        {"SIH3": 4, "SIHO": 1, "SIH2O": 1, "SIH": 1},
        "smiles",
    ),
    (
        "[SiH3][SiH]([SiH3])O[SiH2]O[SiH](O[SiH2]O[SiH]([SiH3])[SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]",  # noqa
        {"SIH3": 6, "SIH": 1, "SIH2O": 3, "SIHO": 3},
        "smiles",
    ),
    (
        "[SiH3][SiH]([SiH3])O[SiH2]O[SiH2][SiH](O[SiH2]O[SiH]([SiH3])[SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]",  # noqa
        {"SIH3": 6, "SIH": 2, "SIH2O": 4, "SIHO": 2},
        "smiles",
    ),
    (
        "[SiH3][SiH]([SiH3])O[SiH2]O[SiH]([SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]",
        {"SIH3": 5, "SIH": 1, "SIH2O": 2, "SIHO": 2},
        "smiles",
    ),
    ("[SiH3]O[SiH2]O[SiH]([SiH3])O[SiH2]O[SiH3]", {}, "smiles"),
    (
        "[SiH3][SiH](O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 7, "SI": 1, "SIHO": 1, "SIO": 1},
        "smiles",
    ),
    (
        "[SiH3][SiH](O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])(O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3])O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3]",  # noqa
        {"SIH3": 13, "SI": 1, "SIHO": 3, "SIO": 3},
        "smiles",
    ),
    (
        "[SiH3][SiH](O[SiH2][Si]([SiH3])(O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3])O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])([SiH3])[SiH3]",  # noqa
        {"SIH3": 13, "SI": 2, "SIH2O": 1, "SIHO": 3, "SIO": 2},
        "smiles",
    ),
    (
        "[SiH3][SiH](O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])([SiH3])O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3]",  # noqa
        {"SIH3": 10, "SI": 1, "SIHO": 2, "SIO": 2},
        "smiles",
    ),
    # Ether concatenation (mix C and Si)
    # triplets
    ("[SiH3]OCO[SiH3]", {}, "smiles"),
    (
        "COCO[SiH2][SiH3]",
        {"CH3O": 1, "CH2O": 1, "SIH2": 1, "SIH3": 1},
        "smiles",
    ),
    ("[SiH3]OCO[SiH2][SiH3]", {"SIH3": 2, "SIH2O": 1, "CH2O": 1}, "smiles"),
    (
        "C[SiH2]OCO[SiH2][SiH3]",
        {"SIH3": 1, "SIH2O": 1, "CH2O": 1, "SIH2": 1, "CH3": 1},
        "smiles",
    ),
    (
        "C[SiH](C)OCO[SiH2][SiH3]",
        {"SIH3": 1, "SIH2O": 1, "CH2O": 1, "SIH": 1, "CH3": 2},
        "smiles",
    ),
    (
        "C[Si](C)(C)OCO[SiH2][SiH3]",
        {"SIH3": 1, "SIH2O": 1, "CH2O": 1, "SI": 1, "CH3": 3},
        "smiles",
    ),
    (
        "CC(C)(C)OCO[SiH2][SiH3]",
        {"SIH3": 1, "SIH2O": 1, "CH2O": 1, "C": 1, "CH3": 3},
        "smiles",
    ),
    (
        "CC(C)OCO[SiH2][SiH3]",
        {"SIH3": 1, "SIH2": 1, "CH2O": 1, "CH-O": 1, "CH3": 2},
        "smiles",
    ),
    (
        "[SiH3]OCO[SiH]([SiH3])[SiH3]",
        {"SIH3": 3, "SIHO": 1, "CH2O": 1},
        "smiles",
    ),
    (
        "C[SiH](C)OCO[SiH]([SiH3])[SiH3]",
        {"SIH3": 2, "CH3": 2, "SIHO": 1, "SIH": 1, "CH2O": 1},
        "smiles",
    ),
    (
        "C[Si](C)(C)OCO[SiH]([SiH3])[SiH3]",
        {"SIH3": 2, "CH3": 3, "SIHO": 1, "SI": 1, "CH2O": 1},
        "smiles",
    ),
    (
        "CC(C)(C)OCO[SiH]([SiH3])[SiH3]",
        {"SIH3": 2, "CH3": 3, "SIHO": 1, "C": 1, "CH2O": 1},
        "smiles",
    ),
    (
        "COCO[SiH]([SiH3])[SiH3]",
        {"SIH3": 2, "CH3O": 1, "SIH": 1, "CH2O": 1},
        "smiles",
    ),
    (
        "CC(C)OCO[SiH]([SiH3])[SiH3]",
        {"SIH3": 2, "CH3": 2, "CH-O": 1, "SIH": 1, "CH2O": 1},
        "smiles",
    ),
    (
        "CCOCO[SiH]([SiH3])[SiH3]",
        {"SIH3": 2, "CH3": 1, "SIH": 1, "CH2O": 2},
        "smiles",
    ),
    (
        "[SiH3]OCO[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 4, "SIO": 1, "CH2O": 1},
        "smiles",
    ),
    (
        "C[Si](C)(C)OCO[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 3, "SIO": 1, "CH2O": 1, "SI": 1, "CH3": 3},
        "smiles",
    ),
    (
        "CC(C)(C)OCO[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 3, "SIO": 1, "CH2O": 1, "C": 1, "CH3": 3},
        "smiles",
    ),
    (
        "COCO[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 3, "SI": 1, "CH2O": 1, "CH3O": 1},
        "smiles",
    ),
    (
        "CC(C)OCO[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 3, "SI": 1, "CH2O": 1, "CH-O": 1, "CH3": 2},
        "smiles",
    ),
    (
        "CC(O[SiH3])O[SiH2][SiH3]",
        {"SIH3": 2, "CH3": 1, "CH-O": 1, "SIH2O": 1},
        "smiles",
    ),
    (
        "CC(O[SiH2][SiH3])O[SiH](C)C",
        {"SIH3": 1, "SIH": 1, "CH3": 3, "CH-O": 1, "SIH2O": 1},
        "smiles",
    ),
    (
        "CC(O[SiH2][SiH3])OC(C)(C)C",
        {"SIH3": 1, "C": 1, "CH3": 4, "CH-O": 1, "SIH2O": 1},
        "smiles",
    ),
    (
        "COC(C)O[SiH2][SiH3]",
        {"SIH3": 1, "CH3": 1, "CH3O": 1, "CH-O": 1, "SIH2": 1},
        "smiles",
    ),
    (
        "CC(C)OC(C)O[SiH2][SiH3]",
        {"SIH3": 1, "CH3": 3, "CH-O": 2, "SIH2": 1},
        "smiles",
    ),
    (
        "CCOC(C)O[SiH2][SiH3]",
        {"SIH3": 1, "CH3": 2, "CH-O": 1, "CH2O": 1, "SIH2": 1},
        "smiles",
    ),
    (
        "CC(O[SiH](C)C)O[SiH]([SiH3])[SiH3]",
        {"SIH3": 2, "CH3": 3, "CH-O": 1, "SIHO": 1, "SIH": 1},
        "smiles",
    ),
    (
        "CC(O[SiH3])O[SiH]([SiH3])[SiH3]",
        {"SIH3": 3, "CH3": 1, "CH-O": 1, "SIHO": 1},
        "smiles",
    ),
    (
        "CC(O[SiH]([SiH3])[SiH3])OC(C)(C)C",
        {"SIH3": 2, "CH3": 4, "CH-O": 1, "SIHO": 1, "C": 1},
        "smiles",
    ),
    (
        "CC(C)OC(C)O[SiH]([SiH3])[SiH3]",
        {"SIH3": 2, "CH3": 3, "CH-O": 2, "SIH": 1},
        "smiles",
    ),
    (
        "COC(C)O[SiH]([SiH3])[SiH3]",
        {"SIH3": 2, "CH3": 1, "CH3O": 1, "CH-O": 1, "SIH": 1},
        "smiles",
    ),
    (
        "CC(O[SiH3])O[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 4, "CH3": 1, "CH-O": 1, "SIO": 1},
        "smiles",
    ),
    (
        "CC(OC(C)(C)C)O[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 3, "CH3": 4, "CH-O": 1, "SIO": 1, "C": 1},
        "smiles",
    ),
    (
        "CC(C)OC(C)O[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 3, "CH3": 3, "CH-O": 2, "SI": 1},
        "smiles",
    ),
    (
        "COC(C)O[Si]([SiH3])([SiH3])[SiH3]",
        {"SIH3": 3, "CH3": 1, "CH3O": 1, "CH-O": 1, "SI": 1},
        "smiles",
    ),
    # quadtruplets
    # TODO: test quadruplets structures, not a priority
    # [SiH3][SiH2]OCOCO[SiH2][SiH3]
]


@pytest.mark.SIO
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sio_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result
