import pytest

from ugropy import get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# Epoxy
# =============================================================================
# PSRK
# =============================================================================
trials_psrk = [
    # propyleneoxide
    ("CC1CO1", {"H2COCH": 1, "CH3": 1}, "smiles"),
    ("C1OC1C1=CC=CC=C1", {"ACH": 5, "AC": 1, "H2COCH": 1}, "smiles"),
    # 2,3-epoxybutane
    ("CC1C(O1)C", {"CH3": 2, "HCOCH": 1}, "smiles"),
    # 2-methyl-2,3-epoxybutane
    ("CC1OC1(C)C", {"CH3": 3, "HCOC": 1}, "smiles"),
    # 2-methyl-1,2-epoxypropane
    ("CC1(CO1)C", {"CH3": 2, "H2COC": 1}, "smiles"),
    # 2,3-dimethyl-2,3-epoxybutane
    ("CC1(C(O1)(C)C)C", {"CH3": 4, "COC": 1}, "smiles"),
]


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_psrk)
def test_51_epoxy_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}


# =============================================================================
# UNIFAC
# =============================================================================
trials_unifac = [
    # propyleneoxide
    ("CC1CO1", {"CH2O": 1, "CH": 1, "CH3": 1}, "smiles"),
    # 2,3-epoxybutane
    ("CC1C(O1)C", {"CH3": 2, "CHO": 1, "CH": 1}, "smiles"),
    # 2-methyl-2,3-epoxybutane
    ("CC1OC1(C)C", {"CH3": 3, "CHO": 1, "C": 1}, "smiles"),
    # 2-methyl-1,2-epoxypropane
    ("CC1(CO1)C", {"CH3": 2, "CH2O": 1, "C": 1}, "smiles"),
    # 2,3-dimethyl-2,3-epoxybutane
    ("CC1(C(O1)(C)C)C", {}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_51_epoxy_unfiac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result

    if mol.subgroups != {}:
        assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}
