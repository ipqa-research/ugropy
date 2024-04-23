import pytest

from ugropy import constantinou_gani_primary, get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# 19- CCN Main group: CH3CN, CH2CN
# =============================================================================

# UNIFAC
trials_unifac = [
    # acetonitrile
    ("CC#N", {"CH3CN": 1}, "smiles"),
    # propionitrile
    ("CCC#N", {"CH3": 1, "CH2CN": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccn_unifac(identifier, result, identifier_type):
    mol = get_groups(unifac, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccn_psrk(identifier, result, identifier_type):
    mol = get_groups(psrk, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ccn_cg(identifier, result, identifier_type):
    mol = get_groups(constantinou_gani_primary, identifier, identifier_type)

    if identifier != "CC#N":
        assert mol.subgroups == result
        assert (
            fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary)
            != {}
        )
    else:
        assert mol.subgroups == {}
