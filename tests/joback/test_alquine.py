import pytest

from ugropy import get_groups, joback
from ugropy.core import fit_atoms


# =============================================================================
# CH, C
# =============================================================================
# Joback
trials = [
    # 1-hexyne
    ("CCCCC#C", {"CH": 1, "C": 1, "-CH2-": 3, "-CH3": 1}, "smiles"),
    (
        "CC#CC1=CC=CC=C1",
        {"ring=CH-": 5, "ring=C<": 1, "C": 2, "-CH3": 1},
        "smiles",
    ),
    # 2-hexyne
    ("CCCC#CC", {"-CH3": 2, "-CH2-": 2, "C": 2}, "smiles"),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_alcohols(identifier, result, identifier_type):
    mol = get_groups(joback, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, joback) != {}
