# =============================================================================
# Description: This file contains the test cases for the particulars molecules
#
# Example, water, pyridine, etc. All the thing that are simple and basic but
# you cannot put anywhere else.
# =============================================================================
from .case import Case


particulars_cases = [
    Case(
        "O",
        "smiles",
        "particulars",
        "Water",
        r=None,
        q=None,
        unifac_result={"H2O": 1},
        psrk_result={"H2O": 1},
        joback_result={},
    ),
]
