import pytest

import ugropy as ug


# =============================================================================
# 31- DOH Main group: DOH
# =============================================================================

# UNIFAC
trials_unifac = [
    # 1,2-ethanediol
    ("C(CO)O", {"DOH": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_doh_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result
