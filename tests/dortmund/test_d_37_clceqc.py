import pytest

import ugropy as ug


# =============================================================================
# 37|CLCC|[69]CL-(C=C)
# =============================================================================

# Dortmund
trials = [
    (
        "ClC(I)=C(Br)C=CC=C",
        {"CH2=CH": 1, "CH=CH": 1, "C=C": 1, "I": 1, "BR": 1, "CL-(C=C)": 1},
        "smiles",
    ),
    # trichloroethylene
    ("C(=C(Cl)Cl)Cl", {"CH=C": 1, "CL-(C=C)": 3}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_cleqc_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result
