import pytest

import ugropy as ug


# =============================================================================
# 35- DMSO Main group: DMSO
# =============================================================================

# UNIFAC
trials_unifac = [
    # dimethyl sulfoxide
    ("CS(=O)C", {"DMSO": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_alquine_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result
