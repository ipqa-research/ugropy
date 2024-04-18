import pytest

from ugropy import constantinou_gani_primary, get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 18- Pyridine Main group: C5H5N, C5H4N, C5H3N
# =============================================================================

# UNIFAC
trials_unifac = [
    # pyridine sharing side with benzene ring
    (
        "C1=CC2=C(C=C1)C1=C(C=NC=C1)C1=C2C=CC=C1",
        {"C5H3N": 1, "ACH": 8, "AC": 4},
        "smiles",
    ),
    # cyclic ester with pyridine
    ("O=C1CCC2=C(O1)C=CN=C2", {"C5H3N": 1, "CH2COO": 1, "CH2": 1}, "smiles"),
    # Pyrydine sharin CH2 with phenyl
    (
        "CCCC1=CC=C(CC2=CC=NC=C2)C=C1",
        {"C5H4N": 1, "ACCH2": 2, "ACH": 4, "CH2": 1, "CH3": 1},
        "smiles",
    ),
    # Nicotine
    (
        "CN1CCCC1C2=CN=CC=C2",
        {"C5H4N": 1, "CH3N": 1, "CH2": 3, "CH": 1},
        "smiles",
    ),
    # Impossible
    ("CC1=CC(C)=C(C)C=N1", {}, "smiles"),
    ("C1=CC=C(C=C1)C1=CC=NC=C1", {"C5H4N": 1, "AC": 1, "ACH": 5}, "smiles"),
    ("CC(=C)C1=CC=NC=C1", {"C5H4N": 1, "CH2=C": 1, "CH3": 1}, "smiles"),
    ("CC1=NC=CC(O)=C1", {"C5H3N": 1, "CH3": 1, "OH": 1}, "smiles"),
    # pyridine
    ("C1=CC=NC=C1", {"C5H5N": 1}, "smiles"),
    # 3-methylpyridine
    ("CC1=CN=CC=C1", {"C5H4N": 1, "CH3": 1}, "smiles"),
    # 2,3-Dimethylpyridine
    ("CC1=C(N=CC=C1)C", {"C5H3N": 1, "CH3": 2}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_cg(identifier, result, identifier_type):
    mol = get_groups(constantinou_gani_primary, identifier, identifier_type)

    if identifier != "C1=CC=NC=C1":
        
        assert mol.subgroups == result

        if mol.subgroups != {}:
            assert fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary) != {}
    else:
        assert mol.subgroups == {}
