import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms

# =============================================================================
# 44- NMP Main group: NMP
# =============================================================================

# UNIFAC
trials_unifac = [
    # N-Methyl-2-pyrrolidone
    ("CN1CCCC1=O", {"NMP": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_npm_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_npm_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}
