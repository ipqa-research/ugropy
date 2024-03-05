import pytest

from ugropy import get_groups, joback
from ugropy.core import fit_atoms


# =============================================================================
# -CH3, -CH2-, >CH-, >C<, CH2=CH-, -CH=CH-, =C<, =C=, CH, C
# =============================================================================
# Joback
trials = [
    (
        "CC=CC(C)=C(C)C=C",
        {"-CH3": 3, "=CH2": 1, "=C<": 2, "=CH-": 3},
        "smiles",
    ),
    (
        "CCCCC=C",
        {"=CH2": 1, "=CH-": 1, "-CH2-": 3, "-CH3": 1},
        "smiles",
    ),  # 1-hexene
    (
        "CCC(CC)C(C)(C)C",
        {"-CH3": 5, "-CH2-": 2, ">CH-": 1, ">C<": 1},
        "smiles",
    ),
    ("CC", {"-CH3": 2}, "smiles"),  # ethane
    ("CCCCCC", {"-CH3": 2, "-CH2-": 4}, "smiles"),  # hexane
    ("CC(C)C", {"-CH3": 3, ">CH-": 1}, "smiles"),  # 2-methylpropane
    ("CC(C)(C)C", {"-CH3": 4, ">C<": 1}, "smiles"),  # 2,2-dimethylpropane
    ("CC(=C(C)C)C", {"-CH3": 4, "=C<": 2}, "smiles"),  # 2,3-dimethylbutene-2
    (
        "CC=C(C)C",
        {"-CH3": 3, "=CH-": 1, "=C<": 1},
        "smiles",
    ),  # 2-methyl-2-butene
    # 2-methyl-1-butene
    ("CCC(=C)C", {"-CH3": 2, "=CH2": 1, "=C<": 1, "-CH2-": 1}, "smiles"),
    ("CCCC=CC", {"-CH2-": 2, "=CH-": 2, "-CH3": 2}, "smiles"),  # 2-hexene
    ("C=C=C", {"=CH2": 2, "=C=": 1}, "smiles"),
    (
        "CC=CC(C)C(C)=C=C",
        {"-CH3": 3, "=CH2": 1, "=C=": 1, "=C<": 1, ">CH-": 1, "=CH-": 2},
        "smiles",
    ),
    ("CC(C)=C=C(C)C", {"-CH3": 4, "=C<": 2, "=C=": 1}, "smiles"),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_no_cyclic_hydrocarbon(identifier, result, identifier_type):
    mol = get_groups(joback, identifier, identifier_type)
    assert mol.subgroups == result
    assert fit_atoms(mol.mol_object, mol.subgroups, joback) != {}
