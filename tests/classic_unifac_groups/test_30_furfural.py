import pytest

import ugropy as ug


# =============================================================================
# 30- furfural Main group: furfural
# =============================================================================

# UNIFAC
trials_unifac = [
    # furfural
    ("C1=COC(=C1)C=O", {"FURFURAL": 1}, "smiles"),
]


@pytest.mark.FURFURAL
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_furfural_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result
