import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 1- CH2 Main group: CH3, CH2, CH, C
# =============================================================================
# UNIFAC
trials_unifac = [
    ("C1CC2CCCC3CCCC1C23", {"CH2": 8, "CH": 4}, "smiles"),
    ("C1C2CCCCC2C2CCCCC12", {"CH2": 9, "CH": 4}, "smiles"),
    ("C1C2CC1CCCC2", {"CH2": 6, "CH": 2}, "smiles"),
    ("C1CCCCCCCC1", {"CH2": 9}, "smiles"),
    ("C1CCCCCCCC1", {"CH2": 9}, "smiles"),
    ("C1C2CC3CC1CC(C2)C3", {"CH2": 6, "CH": 4}, "smiles"),
    ("C12C3C1C1C2C31", {"CH": 6}, "smiles"),
    ("C1CC2CC1CCC2", {"CH2": 6, "CH": 2}, "smiles"),
    ("C1CC2CC3CCC2CC13", {"CH2": 6, "CH": 4}, "smiles"),
    ("C12C3C4C1C1C2C3C41", {"CH": 8}, "smiles"),
    ("C1CC1", {"CH2": 3}, "smiles"),
    ("C1CCC1", {"CH2": 4}, "smiles"),
    ("C1C2CC1CCCC2", {"CH2": 6, "CH": 2}, "smiles"),
    ("CC12C3CCC4CCC1C234", {"CH3": 1, "CH2": 4, "CH": 3, "C": 2}, "smiles"),
    ("CCC(CC)C(C)(C)C", {"CH3": 5, "CH2": 2, "CH": 1, "C": 1}, "smiles"),
    ("C1CCC2CCCCC2C1", {"CH2": 8, "CH": 2}, "smiles"),
    ("C1CCC(CC1)CC2CCCCC2", {"CH2": 11, "CH": 2}, "smiles"),
    ("CC", {"CH3": 2}, "smiles"),  # ethane
    ("CCCCCC", {"CH3": 2, "CH2": 4}, "smiles"),  # hexane
    ("CC(C)C", {"CH3": 3, "CH": 1}, "smiles"),  # 2-methylpropane
    ("CC(C)(C)C", {"CH3": 4, "C": 1}, "smiles"),  # 2,2-dimethylpropane
    ("C1CCCCC1", {"CH2": 6}, "smiles"),  # cyclohexane
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch2_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}
