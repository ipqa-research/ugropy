import pytest

from ugropy import get_groups, psrk, unifac


# =============================================================================
# 16- (C)3N Main group: CH3N, CH2N
# =============================================================================

# UNIFAC
trials_unifac = [
    # Quinuclidine
    ("C1CN2CCC1CC2", {"CH2N": 1, "CH2": 5, "CH": 1}, "smiles"),
    (
        "CCN(CC(=O)CC)C1=CC=CC=C1",
        {"CH3": 2, "CH2N": 1, "AC": 1, "ACH": 5, "CH2CO": 1, "CH2": 1},
        "smiles",
    ),
    ("CCN(C(C)C)C(C)C", {"CH2N": 1, "CH3": 5, "CH": 2}, "smiles"),
    ("CCN(C)CC", {"CH3N": 1, "CH2": 2, "CH3": 2}, "smiles"),
    # trimethylamine
    ("CN(C)C", {"CH3": 2, "CH3N": 1}, "smiles"),
    # triethylamine
    ("CCN(CC)CC", {"CH3": 3, "CH2": 2, "CH2N": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_unifac(identifier, result, identifier_type):
    assert get_groups(unifac, identifier, identifier_type).subgroups == result


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_ch3nh2_psrk(identifier, result, identifier_type):
    assert get_groups(psrk, identifier, identifier_type).subgroups == result
