import pytest

import ugropy as ug


# =============================================================================
# 14 - OH (P), OH (S), OH (T)
# =============================================================================
# Dortmund
trials_dortmund = [

]


@pytest.mark.OH
@pytest.mark.DORTMUND
@pytest.mark.parametrize("identifier, result, identifier_type", trials_dortmund)
def test_unifac_ch2(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    
