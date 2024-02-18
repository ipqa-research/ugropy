import pytest

import ugropy as ug


# =============================================================================
# 30|FURFURAL|[61]FURFURAL
# =============================================================================

# Dortmund
trials = [
    # furfural
    ("C1=COC(=C1)C=O", {"FURFURAL": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_furfural_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result
