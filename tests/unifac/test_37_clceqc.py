import ugropy as ug

import pytest


# =============================================================================
# 37- CL-(C=C) Main group: CL-(C=C)
# =============================================================================

# UNIFAC
trials_unifac = [
    ("ClC(I)=C(Br)C=CC=C", {"CH2=CH": 1, "CH=CH": 1, "C=C": 1, "I": 1, "BR": 1, "CL-(C=C)": 1}, "smiles"),
    ("trichloroethylene", {"CH=C": 1, "CL-(C=C)": 3}, "name"),
]

@pytest.mark.CLCeqC
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cleqc_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result