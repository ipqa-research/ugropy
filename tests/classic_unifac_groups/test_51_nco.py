import pytest

import ugropy as ug


# =============================================================================
# 51- NCO Main group: NCO
# =============================================================================

# UNIFAC
trials_unifac = [
    # Isophorone diisocyanate
    (
        "CC1(CC(CC(C1)(C)CN=C=O)N=C=O)C",
        {"CH3": 3, "CH2": 4, "CH": 1, "C": 2, "NCO": 2},
        "smiles",
    ),
]


@pytest.mark.PSRK
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_thiophene_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == {}
