import pytest

from ugropy import get_groups, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 84- IMIDAZOL Main group: IMIDAZOL
# =============================================================================

# UNIFAC
trials_unifac = [
    # Imidazol
    ("N1C=CN=C1", {"IMIDAZOL": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_thiophene_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}
