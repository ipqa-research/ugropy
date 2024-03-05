import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 26- CNO2 Main group: CH3NO2, CH2NO2, CHNO2
# =============================================================================

# UNIFAC
trials_unifac = [
    (
        "CCCC1=CC=C(C[N+]([O-])=O)C=C1",
        {"ACH": 4, "AC": 1, "CH2NO2": 1, "ACCH2": 1, "CH3": 1, "CH2": 1},
        "smiles",
    ),
    ("[O-][N+](=O)CC1=CC=CC=C1", {"ACH": 5, "AC": 1, "CH2NO2": 1}, "smiles"),
    # nitromethane
    ("C[N+](=O)[O-]", {"CH3NO2": 1}, "smiles"),
    # 1-nitropropane
    ("CCC[N+](=O)[O-]", {"CH3": 1, "CH2": 1, "CH2NO2": 1}, "smiles"),
    # 2-nitropropane
    ("CC(C)[N+](=O)[O-]", {"CH3": 2, "CHNO2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cno2_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cno2_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}
