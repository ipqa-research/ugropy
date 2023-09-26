import pytest

import ugropy as ug


# =============================================================================
#
# =============================================================================
trials_psrk = [
    # 2-propanethiol
    ("CC(C)S", {"CH3": 2, "CHSH": 1}, "smiles"),
    ("CC(S)C1=CC=CC=C1", {"CH3": 1, "CHSH": 1, "AC": 1, "ACH": 5}, "smiles"),
    # 2-methyl-2-propanethiol
    ("CC(C)(C)S", {"CH3": 3, "CSH": 1}, "smiles"),
]


@pytest.mark.PSRK
@pytest.mark.CH3SH
@pytest.mark.parametrize("identifier, result, identifier_type", trials_psrk)
def test_29_ch3sh(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == {}
    assert groups.psrk_groups == result
