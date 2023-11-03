import pytest

import ugropy as ug


# =============================================================================
# 
# =============================================================================
# Joback
trials_unifac = [
    # ("C1CC2CCCC3CCCC1C23", {"CH2": 8, "CH": 4}, "smiles"),
    # ("C1C2CCCCC2C2CCCCC12", {"CH2": 9, "CH": 4}, "smiles"),
    # ("C1C2CC1CCCC2", {"CH2": 6, "CH": 2}, "smiles"),
    # ("C1CCCCCCCC1", {"CH2": 9}, "smiles"),
    # ("C1CCCCCCCC1", {"CH2": 9}, "smiles"),
    # ("C1C2CC3CC1CC(C2)C3", {"CH2": 6, "CH": 4}, "smiles"),
    # ("C12C3C1C1C2C31", {"CH": 6}, "smiles"),
    # ("C1CC2CC1CCC2", {"CH2": 6, "CH": 2}, "smiles"),
    # ("C1CC2CC3CCC2CC13", {"CH2": 6, "CH": 4}, "smiles"),
    # ("C12C3C4C1C1C2C3C41", {"CH": 8}, "smiles"),
    # ("C1CC1", {"CH2": 3}, "smiles"),
    # ("C1CCC1", {"CH2": 4}, "smiles"),
    # ("C1C2CC1CCCC2", {"CH2": 6, "CH": 2}, "smiles"),
    # ("CC12C3CCC4CCC1C234", {"CH3": 1, "CH2": 4, "CH": 3, "C": 2}, "smiles"),
    # ("C1CCC2CCCCC2C1", {"CH2": 8, "CH": 2}, "smiles"),
    # ("C1CCC(CC1)CC2CCCCC2", {"CH2": 11, "CH": 2}, "smiles"),
    # ("C1CCCCC1", {"CH2": 6}, "smiles"),  # cyclohexane
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_joback_cyclic_hydrocarbon(identifier, result, identifier_type):
    assert ug.get_joback_groups(identifier, identifier_type) == result
