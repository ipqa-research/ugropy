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


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_furfural_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result
