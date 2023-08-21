import ugropy as ug

import pytest


# =============================================================================
# A mix of a lot of subgroups
# =============================================================================

# UNIFAC with name
trials_unifac = [

]

@pytest.mark.mix
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_name_unifac(name, result):
    groups = ug.Groups(name)
    assert groups.unifac_groups == result

# UNIFAC with smiles
trials_unifac = [
    # Composed structure problem
    ("C1(=CC=CC=C1)COC(C)(C)C", {"ACH": 5, "AC": 1, "CH2O": 1, "CH3": 3, "C": 1}),
    ("C1(=CC=CC=C1)C(OC(C)(C)C)C", {"ACH": 5, "AC": 1, "CH-O": 1, "CH3": 4, "C": 1}),
    ("C12=CC=CC=C1COC2", {"ACH": 4, "AC": 1, "CH2O": 1, "ACCH2": 1})
]

@pytest.mark.mix
@pytest.mark.UNIFAC
@pytest.mark.parametrize("smiles, result", trials_unifac)
def test_smiles_unifac(smiles, result):
    groups = ug.Groups(smiles, identifier_type="smiles")
    assert groups.unifac_groups == result
