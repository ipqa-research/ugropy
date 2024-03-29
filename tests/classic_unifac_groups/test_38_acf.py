import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 38- ACF Main group: ACF
# =============================================================================

# UNIFAC
trials_unifac = [
    # hexafluorobenzene
    ("C1(=C(C(=C(C(=C1F)F)F)F)F)F", {"ACF": 6}, "smiles"),
    ("FC1=CC=NC=C1", {}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acf_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acf_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}
