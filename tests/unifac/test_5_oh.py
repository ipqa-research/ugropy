import ugropy as ug

import pytest


# =============================================================================
# 5- OH Main group: OH, CH3OH
# =============================================================================

# UNIFAC
trials_unifac = [
    ("1,2-Cyclohexanediol, 4-tert-butyl-1-phenyl-, stereoisomer", {"CH3": 3, "CH2": 3, "CH": 2, "C": 2, "OH": 2, "AC": 1, "ACH": 5}),
    ("(2S,3S)-2-Methyl-1,3-hexanediol", {"CH3": 2, "CH2": 3, "CH": 2, "OH": 2}),
    ("2-propanol", {"CH3": 2, "CH": 1, "OH": 1}),
    ("methanol", {"CH3OH": 1})
]

@pytest.mark.OH
@pytest.mark.UNIFAC
@pytest.mark.parametrize("name, result", trials_unifac)
def test_ch2_unifac(name, result):
    substance = ug.Substance(name)
    assert substance.unifac_groups == result