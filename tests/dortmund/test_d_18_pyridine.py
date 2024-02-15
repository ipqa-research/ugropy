import pytest

import ugropy as ug


# =============================================================================
# 18|PYRIDINE|[37]AC2H2N [38]AC2HN [39]AC2N
# =============================================================================

# Dortmund
trials = [
    ("CC1=CC=CC(O)=N1", {"AC2N": 1, "ACH": 3, "CH3": 1, "OH": 1}, "smiles"),
    ("CC1=CC=CC(C)=N1", {"AC2N": 1, "ACH": 3, "CH3": 2}, "smiles"),
    # pyridine sharing side with benzene ring
    (
       "C1=CC2=C(C=C1)C1=C(C=NC=C1)C1=C2C=CC=C1",
       {"ACH": 9, "AC": 6, "AC2H2N": 1},
       "smiles",
    ),
    # Pyrydine sharin CH2 with phenyl
    (
       "CCCC1=CC=C(CC2=CC=NC=C2)C=C1",
       {"ACH": 6, "AC2H2N": 1, "ACCH2": 2, "AC": 1, "CH2": 1, "CH3": 1},
       "smiles",
    ),
    ("CC1=CC(C)=C(C)C=N1", {"ACH": 1, "AC2HN": 1, "ACCH3": 2, "CH3": 1}, "smiles"),
    ("C1=CC=C(C=C1)C1=CC=NC=C1", {"AC2H2N": 1, "AC": 2, "ACH": 7}, "smiles"),
    ("CC(=C)C1=CC=NC=C1", {"ACH": 2, "AC2H2N": 1, "AC": 1, "CH2=C": 1, "CH3": 1}, "smiles"),
    ("CC1=NC=CC(O)=C1", {"ACH": 2, "AC2HN": 1, "ACOH": 1, "CH3": 1}, "smiles"),
    # pyridine
    ("C1=CC=NC=C1", {"ACH": 3, "AC2H2N": 1}, "smiles"),
    # 3-methylpyridine
    ("CC1=CN=CC=C1", {"ACH": 2, "AC2H2N": 1, "ACCH3": 1}, "smiles"),
    # 2,3-Dimethylpyridine
    ("CC1=C(N=CC=C1)C", {"ACH": 2, "AC2HN": 1, "ACCH3": 1, "CH3": 1}, "smiles"),
]


@pytest.mark.Dortmund
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_ch3nh2_dortmund(identifier, result, identifier_type):
    assert ug.get_dortmund_groups(identifier, identifier_type) == result
