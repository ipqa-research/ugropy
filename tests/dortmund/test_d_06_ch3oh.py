import pytest

import ugropy as ug


# =============================================================================
# 6|CH3OH|[15]CH3OH
# =============================================================================

# Dortmund
# methanol
trials = [("CO", {"CH3OH": 1}, "smiles")]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_oh_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result
