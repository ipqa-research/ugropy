import pytest

import ugropy as ug


# =============================================================================
# 39- DMF Main group: DMF, HCON(CH2)2
# =============================================================================

# UNIFAC
trials_unifac = [
    # N,N-Dimethylformamide
    ("CN(C)C=O", {"DMF": 1}, "smiles"),
    (
        "BrCN(CC1=CC=NC=C1)C=O",
        {"HCON(CH2)2": 1, "BR": 1, "C5H4N": 1},
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


@pytest.mark.DMF
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_dmf_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result
