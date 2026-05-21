import pytest

import ugropy as ug

# Smiles for testing
eter = "CCOC(C)C"
ester = "COC(=O)CC1=CC=CC=C1"


def test_filter_mostly_polarity():
    sol1 = ug.unifac.get_groups(
        ester, "smiles", search_multiple_solutions=True
    )
    fil1 = ug.unifac.filter_mostly_polarity(sol1)

    assert fil1[0].subgroups == {"CH3": 1, "ACH": 5, "AC": 1, "CH2COO": 1}

    fil2 = ug.unifac.filter_mostly_polarity(sol1, "apolar")

    assert fil2[0].subgroups == {"CH3": 1, "ACH": 5, "ACCH2": 1, "COO": 1}

    with pytest.raises(ValueError):
        ug.unifac.filter_mostly_polarity(sol1, "invalid")


def test_filter_polyatomic_criteria():
    sol1 = ug.unifac.get_groups(eter, "smiles", search_multiple_solutions=True)
    fil1 = ug.unifac.filter_polyatomic_criteria(sol1)

    assert fil1[0].subgroups == {"CH3": 3, "CH2O": 1, "CH": 1}

    sol2 = ug.unifac.get_groups(
        ester, "smiles", search_multiple_solutions=True
    )
    fil2 = ug.unifac.filter_polyatomic_criteria(sol2)

    assert fil2[0].subgroups == {"CH3": 1, "ACH": 5, "ACCH2": 1, "COO": 1}

    with pytest.raises(ValueError):
        ug.unifac.filter_polyatomic_criteria(sol1, "F")


def test_filter_polarity_contribution():
    sol1 = ug.unifac.get_groups(
        ester, "smiles", search_multiple_solutions=True
    )
    fil1 = ug.unifac.filter_polarity_contribution(sol1, "Q", "polar")

    assert fil1[0].subgroups == {"CH3": 1, "ACH": 5, "AC": 1, "CH2COO": 1}

    fil2 = ug.unifac.filter_polarity_contribution(sol1, "Q", "apolar")

    assert fil2[0].subgroups == {"CH3": 1, "ACH": 5, "ACCH2": 1, "COO": 1}

    with pytest.raises(ValueError):
        ug.unifac.filter_polarity_contribution(sol1, "Q", "invalid")

    with pytest.raises(ValueError):
        ug.unifac.filter_polarity_contribution(sol1, "F", "polar")
