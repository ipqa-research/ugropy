import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 20- COOH Main group: COOH, HCOOH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("OC(O)=O", {"COOH": 1, "OH": 1}, "smiles"),
    ("CCOC(O)=O", {"COOH": 1, "CH2O": 1, "CH3": 1}, "smiles"),
    # 2,4-Diaminobutyric acid
    (
        "C(CN)C(C(=O)O)N",
        {"COOH": 1, "CHNH2": 1, "CH2": 1, "CH2NH2": 1},
        "smiles",
    ),
    # acetic acid
    ("CC(=O)O", {"CH3": 1, "COOH": 1}, "smiles"),
    # formic acid
    ("C(=O)O", {"HCOOH": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cooh_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cooh_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}
