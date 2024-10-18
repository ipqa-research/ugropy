# =============================================================================
# Description: This file contains the test cases for the acids molecules.
# =============================================================================
from .case import Case


acids_cases = [
    Case(
        identifier="OC(O)=O",
        identifier_type="smiles",
        cases_module="acids",
        unifac_result={"OH": 1, "COOH": 1},
        psrk_result={"OH": 1, "COOH": 1},
        joback_result={"-OH (alcohol)": 1, "-COOH (acid)": 1},
    ),
    Case(
        identifier="CCOC(O)=O",
        identifier_type="smiles",
        cases_module="acids",
        unifac_result={"CH3": 1, "CH2O": 1, "COOH": 1},
        psrk_result={"CH3": 1, "CH2O": 1, "COOH": 1},
        joback_result=[
            {"-CH3": 1, "-CH2-": 1, "-O- (non-ring)": 1, "-COOH (acid)": 1},
            {"-CH3": 1, "-CH2-": 1, "-OH (alcohol)": 1, "-COO- (ester)": 1},
        ],
    ),
    Case(
        identifier="C(CN)C(C(=O)O)N",
        identifier_type="smiles",
        cases_module="acids",
        unifac_result={"CH2": 1, "CH2NH2": 1, "CHNH2": 1, "COOH": 1},
        psrk_result={"CH2": 1, "CH2NH2": 1, "CHNH2": 1, "COOH": 1},
        joback_result={"-CH2-": 2, ">CH-": 1, "-COOH (acid)": 1, "-NH2": 2},
    ),
    Case(
        identifier="CC(=O)O",
        identifier_type="smiles",
        cases_module="acids",
        unifac_result={"CH3": 1, "COOH": 1},
        psrk_result={"CH3": 1, "COOH": 1},
        joback_result={"-CH3": 1, "-COOH (acid)": 1},
    ),
    Case(
        identifier="C(=O)O",
        identifier_type="smiles",
        cases_module="acids",
        unifac_result={"HCOOH": 1},
        psrk_result={"HCOOH": 1},
        joback_result={},
    ),
]
