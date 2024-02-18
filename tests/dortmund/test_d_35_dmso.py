import pytest

import ugropy as ug


# =============================================================================
# 35|DMSO|[67]DMSO
# =============================================================================

# Dortmund
trials = [
    # dimethyl sulfoxide
    ("CS(=O)C", {"DMSO": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_alquine_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result
