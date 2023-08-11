import ugropy as ug

import pytest


# =============================================================================
# 46- CON Main group: AMH2, AMHCH3, AMHCH2, AM(CH3)2, AMCH3CH2, AM(CH2)2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("acetamide", {"CH3": 1, "AMH2": 1}, "name"),
    ("N-Methylacetamide", {"CH3": 1, "AMHCH3": 1}, "name"),
    ("N-Ethylacetamide", {"CH3": 2, "AMHCH2": 1}, "name"),
    ("N,N-Dimethylacetamide", {"CH3": 1, "AM(CH3)2": 1}, "name"),
    ("N-ethyl-N-methylacetamide", {"CH3": 2, "AMCH3CH2": 1}, "name"),
    ("N,N-Diethylacetamide", {"CH3": 3, "AM(CH2)2": 1}, "name"),
]

@pytest.mark.CON
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_con_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result