import pytest

import ugropy as ug


# =============================================================================
# 37- CL-(C=C) Main group: CL-(C=C)
# =============================================================================

# UNIFAC
trials_unifac = [
    (
        "ClC(I)=C(Br)C=CC=C",
        {"CH2=CH": 1, "CH=CH": 1, "C=C": 1, "I": 1, "BR": 1, "CL-(C=C)": 1},
        "smiles",
    ),
    # trichloroethylene
    ("C(=C(Cl)Cl)Cl", {"CH=C": 1, "CL-(C=C)": 3}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cleqc_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result
