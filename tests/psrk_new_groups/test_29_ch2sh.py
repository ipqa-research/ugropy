import pytest

from ugropy import get_groups, psrk
from ugropy.core import fit_atoms


# =============================================================================
#
# =============================================================================
trials_psrk = [
    # 2-propanethiol
    ("CC(C)S", {"CH3": 2, "CHSH": 1}, "smiles"),
    ("CC(S)C1=CC=CC=C1", {"CH3": 1, "CHSH": 1, "AC": 1, "ACH": 5}, "smiles"),
    # 2-methyl-2-propanethiol
    ("CC(C)(C)S", {"CH3": 3, "CSH": 1}, "smiles"),
]


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_psrk)
def test_29_ch3sh(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}
