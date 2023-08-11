import ugropy as ug

import pytest


# =============================================================================
# 44- NPM Main group: NPM
# =============================================================================

# UNIFAC
trials_unifac = [
    ("N-Methyl-2-pyrrolidone", {"NMP": 1}, "name"),
]

@pytest.mark.NMP
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_npm_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result