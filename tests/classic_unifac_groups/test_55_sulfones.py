import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 55- sulfones Main group: (CH2)2SU, CH2CHSU
# =============================================================================

# UNIFAC
trials_unifac = [
    # sulfolane
    ("C1CCS(=O)(=O)C1", {"(CH2)2SU": 1, "CH2": 2}, "smiles"),
    # 2,4-dimethylsulfolane
    (
        "CC1CC(S(=O)(=O)C1)C",
        {"CH2CHSU": 1, "CH3": 2, "CH2": 1, "CH": 1},
        "smiles",
    ),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sulfones_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sulfones_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == {}
