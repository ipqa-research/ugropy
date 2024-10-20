# =============================================================================
# Epoxides test cases
# =============================================================================
from .case import Case


epoxides_cases = [
    Case(
        identifier="O1C2OC12",
        identifier_type="smiles",
        cases_module="epoxides",
        unifac_result={"CHO": 2},
        psrk_result={"CHO": 2},
        joback_result={"-O- (ring)": 2, "ring>CH-": 2},
    ),
    Case(
        identifier="O1C23OC22OC132",
        identifier_type="smiles",
        cases_module="epoxides",
        unifac_result={},
        psrk_result={},
        joback_result={"-O- (ring)": 3, "ring>C<": 3},
    ),
    Case(
        identifier="C1OC11CO1",
        identifier_type="smiles",
        cases_module="epoxides",
        unifac_result={"CH2O": 2, "C": 1},
        psrk_result={"H2COC": 1, "CH2O": 1},
        joback_result={"-O- (ring)": 2, "ring>C<": 1, "ring-CH2-": 2},
    ),
    Case(
        identifier="CC1CO1",
        identifier_type="smiles",
        cases_module="epoxides",
        unifac_result=[
            {"CH3": 1, "CH": 1, "CH2O": 1},
            {"CH3": 1, "CH2": 1, "CHO": 1},
        ],
        psrk_result={"CH3": 1, "H2COCH": 1},
        joback_result={
            "-CH3": 1,
            "ring-CH2-": 1,
            "ring>CH-": 1,
            "-O- (ring)": 1,
        },
    ),
    Case(
        identifier="C1OC1C1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="epoxides",
        unifac_result={"ACH": 5, "ACCH": 1, "CH2O": 1},
        psrk_result=[
            {"ACH": 5, "AC": 1, "H2COCH": 1},
            {"ACH": 5, "ACCH": 1, "CH2O": 1},
        ],
        joback_result={
            "ring-CH2-": 1,
            "ring>CH-": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-O- (ring)": 1,
        },
    ),
    Case(
        identifier="CC1C(O1)C",
        identifier_type="smiles",
        cases_module="epoxides",
        unifac_result={"CH3": 2, "CH": 1, "CHO": 1},
        psrk_result={"CH3": 2, "HCOCH": 1},
        joback_result={"-CH3": 2, "ring>CH-": 2, "-O- (ring)": 1},
    ),
    Case(
        identifier="CC1OC1(C)C",
        identifier_type="smiles",
        cases_module="epoxides",
        unifac_result={"CH3": 3, "C": 1, "CHO": 1},
        psrk_result={"CH3": 3, "HCOC": 1},
        joback_result={
            "-CH3": 3,
            "ring>CH-": 1,
            "ring>C<": 1,
            "-O- (ring)": 1,
        },
    ),
    Case(
        identifier="CC1(CO1)C",
        identifier_type="smiles",
        cases_module="epoxides",
        unifac_result={"CH3": 2, "C": 1, "CH2O": 1},
        psrk_result={"CH3": 2, "H2COC": 1},
        joback_result={
            "-CH3": 2,
            "ring-CH2-": 1,
            "ring>C<": 1,
            "-O- (ring)": 1,
        },
    ),
    Case(
        identifier="CC1(C(O1)(C)C)C",
        identifier_type="smiles",
        cases_module="epoxides",
        unifac_result={},
        psrk_result={"CH3": 4, "COC": 1},
        joback_result={"-CH3": 4, "ring>C<": 2, "-O- (ring)": 1},
    ),
]
