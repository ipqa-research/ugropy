import pytest

import ugropy as ug

# =============================================================================
# 44- NMP Main group: NMP
# =============================================================================

# UNIFAC
trials_unifac = [
    # N-Methyl-2-pyrrolidone
    ("CN1CCCC1=O", {"NMP": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_npm_unifac(identifier, result, identifier_type):
    assert ug.get_groups(ug.unifac, identifier, identifier_type) == result
    assert ug.get_groups(ug.psrk, identifier, identifier_type) == result
