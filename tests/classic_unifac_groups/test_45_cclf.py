import pytest

import ugropy as ug


# =============================================================================
# 45- CCLF Main group: CCL3F, CCL2F, HCCL2F, HCCLF, CCLF2, HCCLF2, CCLF3,
#                      CCL2F2
# =============================================================================

# UNIFAC
trials_unifac = [
    ("Trichlorofluoromethane", {"CCL3F": 1}, "name"),
    ("Tetrachloro-1,2-difluoroethane", {"CCL2F": 2}, "name"),
    ("Dichlorofluoromethane", {"HCCL2F": 1}, "name"),
    ("1-Chloro-1,2,2,2-tetrafluoroethane", {"CF3": 1, "HCCLF": 1}, "name"),
    ("1,2-Dichlorotetrafluoroethane", {"CCLF2": 2}, "name"),
    ("Chlorodifluoromethane", {"HCCLF2": 1}, "name"),
    ("Chlorotrifluoromethane", {"CCLF3": 1}, "name"),
    ("Dichlorodifluoromethane", {"CCL2F2": 1}, "name"),
]


@pytest.mark.CCLF
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cclf_unifac(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type)
    assert groups.unifac_groups == result
