import pytest

import ugropy as ug


# =============================================================================
# -CH3, -CH2-, >CH-, >C<, CH2=CH-, -CH=CH-, =C<, =C=, CH, C
# =============================================================================
# Joback
trials_unifac = [
    ("CCC(CC)C(C)(C)C", {"-CH3": 5, "-CH2-": 2, ">CH-": 1, ">C<": 1}, "smiles"),
    ("CC", {"-CH3": 2}, "smiles"),  # ethane
    ("CCCCCC", {"-CH3": 2, "-CH2-": 4}, "smiles"),  # hexane
    ("CC(C)C", {"-CH3": 3, ">CH-": 1}, "smiles"),  # 2-methylpropane
    ("CC(C)(C)C", {"-CH3": 4, ">C<": 1}, "smiles"),  # 2,2-dimethylpropane
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_joback_no_cyclic_hydrocarbon(identifier, result, identifier_type):
    assert ug.get_joback_groups(identifier, identifier_type) == result
