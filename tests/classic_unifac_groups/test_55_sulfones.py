import pytest

import ugropy as ug


# =============================================================================
# 55- sulfones Main group: (CH2)2SU, CH2CHSU
# =============================================================================

# UNIFAC
trials_unifac = [
    # sulfolane
    ("C1CCS(=O)(=O)C1", {"(CH2)2SU": 1, "CH2": 2}, "smiles"),
    # 2,4-dimethylsulfolane
    (
        "CC1CC(S(=O)(=O)C1)C",
        {"CH2CHSU": 1, "CH3": 2, "CH2": 1, "CH": 1},
        "smiles",
    ),
]


@pytest.mark.sulfones
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_sulfones_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
    assert groups.psrk_groups == {}
