import pytest

import ugropy as ug


# =============================================================================
# 38- ACF Main group: ACF
# =============================================================================

# UNIFAC
trials_unifac = [
    # hexafluorobenzene
    ("C1(=C(C(=C(C(=C1F)F)F)F)F)F", {"ACF": 6}, "smiles"),
    ("FC1=CC=NC=C1", {}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acf_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result
