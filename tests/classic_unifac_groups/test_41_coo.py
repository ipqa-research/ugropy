import pytest

from ugropy import constantinou_gani_primary, get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 41- COO Main group: COO
# =============================================================================

# UNIFAC
trials_unifac = [
    # Ascorbic acid
    (
        "OCC(O)C1OC(=O)C(O)=C1O",
        {"COO": 1, "C=C": 1, "OH": 4, "CH": 2, "CH2": 1},
        "smiles",
    ),
    # Procaine
    (
        "CCN(CC)CCOC(=O)C1=CC=C(N)C=C1",
        {
            "ACNH2": 1,
            "ACH": 4,
            "AC": 1,
            "COO": 1,
            "CH2": 3,
            "CH3": 2,
            "CH2N": 1,
        },
        "smiles",
    ),
    # Cocaine
    (
        "COC(=O)C1C2CCC(CC1OC(=O)C1=CC=CC=C1)N2C",
        {
            "CH3": 1,
            "CH2": 3,
            "CH": 4,
            "CH3N": 1,
            "AC": 1,
            "ACH": 5,
            "COO": 2,
            "CH3N": 1,
        },
        "smiles",
    ),
    # Methyl acrylate
    ("COC(=O)C=C", {"CH3": 1, "CH2=CH": 1, "COO": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_coo_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_coo_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_coo_cg(identifier, result, identifier_type):
    mol = get_groups(constantinou_gani_primary, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary) != {}
