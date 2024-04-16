import pytest

from ugropy import constantinou_gani_primary, get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 2- C=C Main group: CH2=CH, CH=CH, CH2=C, CH=C, C=C
# =============================================================================

# UNIFAC
trials_unifac = [
    (
        "CC=C(C)C1=C(C=CC=C1C=C)C(C)=C(C)C",
        {"CH3": 5, "CH2=CH": 1, "CH=C": 1, "C=C": 1, "ACH": 3, "AC": 3},
        "smiles",
    ),
    (
        "CC=CC(C)=C(C)C=C",
        {"CH2=CH": 1, "CH=CH": 1, "C=C": 1, "CH3": 3},
        "smiles",
    ),
    # alpha-pinene
    (
        "CC1=CCC2CC1C2(C)C",
        {"CH3": 3, "CH2": 2, "CH": 2, "C": 1, "CH=C": 1},
        "smiles",
    ),
    # d-limonene
    (
        "CC1=CCC(CC1)C(=C)C",
        {"CH3": 2, "CH2": 3, "CH": 1, "CH2=C": 1, "CH=C": 1},
        "smiles",
    ),
    # 2,3-Dimethyl-1,3-cyclohexadiene
    ("CC1=CCCC=C1C", {"CH2": 2, "CH=C": 2, "CH3": 2}, "smiles"),
    # 3,3'-(Pentane-1,3-diyl)dicyclohexene
    (
        "CCC(CCC1CCCC=C1)C2CCCC=C2",
        {"CH3": 1, "CH2": 9, "CH": 3, "CH=CH": 2},
        "smiles",
    ),
    ("C1CCC=CC1", {"CH2": 4, "CH=CH": 1}, "smiles"),  # cyclohexene
    ("CC(=C(C)C)C", {"C=C": 1, "CH3": 4}, "smiles"),  # 2,3-dimethylbutene-2
    ("CC=C(C)C", {"CH=C": 1, "CH3": 3}, "smiles"),  # 2-methyl-2-butene
    # 2-methyl-1-butene
    ("CCC(=C)C", {"CH2": 1, "CH2=C": 1, "CH3": 2}, "smiles"),
    ("CCCC=CC", {"CH2": 2, "CH=CH": 1, "CH3": 2}, "smiles"),  # 2-hexene
    ("CCCCC=C", {"CH2": 3, "CH2=CH": 1, "CH3": 1}, "smiles"),  # 1-hexene
    # impossibles
    ("C=C=C", {}, "smiles"),
    ("CC=CC(C)C(C)=C=C", {}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ceqc_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ceqc_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ceqc_cg(identifier, result, identifier_type):
    mol = get_groups(constantinou_gani_primary, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary) != {}
