import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 42- SIH2 Main group: SIH3, SIH2, SIH, SI
# =============================================================================

# UNIFAC
trials_unifac = [
    # methylsilane
    ("C[SiH3]", {"CH3": 1, "SIH3": 1}, "smiles"),
    ("CC[Si](CC)([H])[H]", {"CH3": 2, "CH2": 2, "SIH2": 1}, "smiles"),
    (
        "C[Si](O[Si](C)(C)C)(O[Si](C)(C)C)[H]",
        {"CH3": 7, "SIO": 1, "SIHO": 1, "SI": 1},
        "smiles",
    ),
    # Hexamethyldisiloxane
    ("C[Si](C)(C)O[Si](C)(C)C", {"CH3": 6, "SIO": 1, "SI": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sih2_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sih2_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}
