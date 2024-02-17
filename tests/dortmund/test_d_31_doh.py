import pytest

import ugropy as ug


# =============================================================================
# 31|DOH|[62]DOH
# =============================================================================

# Dortmund
trials = [
    # 1,2-ethanediol
    ("C(CO)O", {"DOH": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_doh_dortmund(identifier, result, identifier_type):
    assert ug.get_groups(ug.dortmund, identifier, identifier_type) == result
