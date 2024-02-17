import pytest

import ugropy as ug


# =============================================================================
# 7|H2O|[16]H2O
# =============================================================================

# Dortmund
trials = [("O", {"H2O": 1}, "smiles")]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_h2o_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result
