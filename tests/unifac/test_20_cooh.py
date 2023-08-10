import ugropy as ug

import pytest


# =============================================================================
# 20- COOH Main group: COOH, HCOOH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("acetic acid", {"CH3": 1, "COOH": 1}, "name"),
    ("formic acid", {"HCOOH": 1}, "name"),
]

@pytest.mark.COOH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cooh_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result