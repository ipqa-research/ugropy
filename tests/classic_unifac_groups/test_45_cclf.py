import pytest

import ugropy as ug


# =============================================================================
# 45- CCLF Main group: CCL3F, CCL2F, HCCL2F, HCCLF, CCLF2, HCCLF2, CCLF3,
#                      CCL2F2
# =============================================================================

# UNIFAC
trials_unifac = [
    # Trichlorofluoromethane
    ("C(F)(Cl)(Cl)Cl", {"CCL3F": 1}, "smiles"),
    # Tetrachloro-1,2-difluoroethane
    ("C(C(F)(Cl)Cl)(F)(Cl)Cl", {"CCL2F": 2}, "smiles"),
    # Dichlorofluoromethane
    ("C(F)(Cl)Cl", {"HCCL2F": 1}, "smiles"),
    # 1-Chloro-1,2,2,2-tetrafluoroethane
    ("C(C(F)(F)F)(F)Cl", {"CF3": 1, "HCCLF": 1}, "smiles"),
    # 1,2-Dichlorotetrafluoroethane
    ("C(C(F)(F)Cl)(F)(F)Cl", {"CCLF2": 2}, "smiles"),
    # Chlorodifluoromethane
    ("C(F)(F)Cl", {"HCCLF2": 1}, "smiles"),
    # Chlorotrifluoromethane
    ("C(F)(F)(F)Cl", {"CCLF3": 1}, "smiles"),
    # Dichlorodifluoromethane
    ("C(F)(F)(Cl)Cl", {"CCL2F2": 1}, "smiles"),
]


@pytest.mark.PSRK
@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_cclf_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
    assert ug.get_psrk_groups(identifier, identifier_type) == result
