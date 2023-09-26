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


@pytest.mark.ACF
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_acf_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result
