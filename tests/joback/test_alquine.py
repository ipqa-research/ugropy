import pytest

import ugropy as ug


# =============================================================================
# CH, C
# =============================================================================
# Joback
trials = [
    (
        "CC#CC1=CC=CC=C1",
        {"ring=CH-": 5, "ring=C<": 1, "C": 2, "-CH3": 1},
        "smiles",
    ),
    # 1-hexyne
    ("CCCCC#C", {"CH": 1, "C": 1, "-CH2-": 3, "-CH3": 1}, "smiles"),
    # 2-hexyne
    ("CCCC#CC", {"-CH3": 2, "-CH2-": 2, "C": 2}, "smiles"),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_alcohols(identifier, result, identifier_type):
    assert ug.get_joback_groups(identifier, identifier_type) == result
