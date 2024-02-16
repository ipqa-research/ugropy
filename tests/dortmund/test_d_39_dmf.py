import pytest

import ugropy as ug


# =============================================================================
# 39|DMF|[72]DMF [73]HCON(CH2)2
# =============================================================================

# Dortmund
trials = [
    # N,N-Dimethylformamide
    ("CN(C)C=O", {"DMF": 1}, "smiles"),
    (
        "BrCN(CC1=CC=NC=C1)C=O",
        {"HCON(CH2)2": 1, "BR": 1, "AC2H2N": 1, "ACH": 2, "AC": 1},
        "smiles",
    ),
    (
        "BrCN(CC1=CC=CC=C1)C=O",
        {"HCON(CH2)2": 1, "BR": 1, "AC": 1, "ACH": 5},
        "smiles",
    ),
    # N,N-Diethylformamide
    ("CCN(CC)C=O", {"CH3": 2, "HCON(CH2)2": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_dmf_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result
