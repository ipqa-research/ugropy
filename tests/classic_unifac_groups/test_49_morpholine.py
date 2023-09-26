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


@pytest.mark.morpholine
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_morpholine_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == result