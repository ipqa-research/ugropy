import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 12- HCOO Main group: HCOO
# =============================================================================

# UNIFAC
trials_unifac = [
    # phenyl formate
    ("C1=CC=C(C=C1)OC=O", {"ACH": 5, "AC": 1, "HCOO": 1}, "smiles"),
    # methyl formate
    ("COC=O", {"HCOO": 1, "CH3": 1}, "smiles"),
    # ethyl formate
    ("CCOC=O", {"HCOO": 1, "CH2": 1, "CH3": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cho_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cho_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}
