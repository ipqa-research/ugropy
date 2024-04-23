import pytest

from ugropy import constantinou_gani_primary, get_groups, psrk, unifac
from ugropy.core import fit_atoms


# =============================================================================
# Composed structures
# =============================================================================

# UNIFAC
trials_unifac = [
    # TODO Thanos explore and test more this kind of monsters
    (
        "CCCC1=C(COC(C)(C)COC(=O)OCC)C=C(CC2=CC=CC=C2)C=C1",
        {
            "CH3": 4,
            "C": 1,
            "ACH": 8,
            "ACCH2": 2,
            "CH2O": 2,
            "COO": 1,
            "CH2": 2,
            "AC": 2,
        },
        "smiles",
    ),
    (
        "CC(C)CC1=CC=C(C=C1)C(C)OC(C)(C)C",
        {"CH3": 6, "CH": 1, "C": 1, "ACH": 4, "ACCH2": 1, "AC": 1, "CHO": 1},
        "smiles",
    ),
    # two solutions
    (
        "CCCC1=CC=C(CC(=O)OC)C=C1",
        [
            {"CH3": 2, "CH2": 1, "ACH": 4, "ACCH2": 1, "AC": 1, "CH2COO": 1},
            {"CH3": 2, "CH2": 1, "ACH": 4, "ACCH2": 2, "COO": 1},
        ],
        "smiles",
    ),
    (
        "C1=CC(=CC=C1COC(C)(C)C)CCC",
        {"ACH": 4, "ACCH2": 1, "AC": 1, "CH2O": 1, "CH3": 4, "CH2": 1, "C": 1},
        "smiles",
    ),
    (
        "C13=C(C=C(C=C1)CC2=CC=CC(=C2)CC)CCCC3",
        {"CH3": 1, "CH2": 2, "ACH": 7, "ACCH2": 4, "AC": 1},
        "smiles",
    ),
    (
        "C13=C(C(=C(C(=C1C)C)CC2=C(C(=C(C(=C2C)CC)O[H])N([H])[H])C)C)CCCC3",
        {
            "CH3": 1,
            "CH2": 2,
            "AC": 1,
            "ACCH3": 5,
            "ACCH2": 4,
            "ACOH": 1,
            "ACNH2": 1,
        },
        "smiles",
    ),
    # toluene - O - tertbutyl (eter between toluene and tertbutyl)
    (
        "C1(=CC=CC=C1)COC(C)(C)C",
        {"ACH": 5, "AC": 1, "CH2O": 1, "C": 1, "CH3": 3},
        "smiles",
    ),
    # diphenyl methane
    ("C1=CC=C(C=C1)CC2=CC=CC=C2", {"ACH": 10, "ACCH2": 1, "AC": 1}, "smiles"),
    (
        "C1(=CC=CC=C1)COC(C)(C)C",
        {"ACH": 5, "AC": 1, "CH2O": 1, "CH3": 3, "C": 1},
        "smiles",
    ),
    (
        "C1(=CC=CC=C1)C(OC(C)(C)C)C",
        {"ACH": 5, "AC": 1, "CHO": 1, "CH3": 4, "C": 1},
        "smiles",
    ),
    ("C12=CC=CC=C1COC2", {"ACH": 4, "AC": 1, "CH2O": 1, "ACCH2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_composed_unifac(identifier, result, identifier_type):
    unifac_result = get_groups(unifac, identifier, identifier_type)
    try:
        mol = get_groups(unifac, identifier, identifier_type)
        assert mol.subgroups == result
        assert fit_atoms(mol.mol_object, mol.subgroups, unifac) != {}
    except ValueError:
        for uni_r, sol in zip(unifac_result.subgroups, result):
            assert uni_r == sol
            assert fit_atoms(mol.mol_object, uni_r, unifac) != {}


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_composed_psrk(identifier, result, identifier_type):
    psrk_result = get_groups(psrk, identifier, identifier_type)
    try:
        mol = get_groups(psrk, identifier, identifier_type)
        assert mol.subgroups == result
        assert fit_atoms(mol.mol_object, mol.subgroups, psrk) != {}
    except ValueError:
        for psrk_r, sol in zip(psrk_result.subgroups, result):
            assert psrk_r == sol
            assert fit_atoms(mol.mol_object, psrk_r, psrk) != {}


@pytest.mark.ConstantinouGani
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_composed_cg(identifier, result, identifier_type):
    cg_result = get_groups(
        constantinou_gani_primary, identifier, identifier_type
    )
    try:
        mol = get_groups(
            constantinou_gani_primary, identifier, identifier_type
        )
        assert mol.subgroups == result
        assert (
            fit_atoms(mol.mol_object, mol.subgroups, constantinou_gani_primary)
            != {}
        )
    except ValueError:
        for psrk_r, sol in zip(cg_result.subgroups, result):
            assert psrk_r == sol
            assert (
                fit_atoms(mol.mol_object, psrk_r, constantinou_gani_primary)
                != {}
            )
