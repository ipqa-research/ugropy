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

]

@pytest.mark.mix
@pytest.mark.UNIFAC
@pytest.mark.parametrize("smiles, result", trials_unifac)
def test_smiles_unifac(smiles, result):
    groups = ug.Groups(smiles, identifier_type="smiles")
    assert groups.unifac_groups == result
