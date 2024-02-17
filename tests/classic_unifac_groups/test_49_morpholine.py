import pytest

import ugropy as ug


# =============================================================================
# 49- morpholine Main group: MORPH
# =============================================================================

# UNIFAC
trials_unifac = [
    # morpholine
    ("C1COCCN1", {"MORPH": 1}, "smiles"),
]


@pytest.mark.PSRK
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_morpholine_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result
