# =============================================================================
# Description: This file contains the test cases for the alcohols module.
#
# Just simple alcohols.
# =============================================================================
from .case import Case


alcohols_cases = [
    Case(
        identifier="C(CO)O",
        identifier_type="smiles",
        cases_module="alcohols",
        r=None,
        q=None,
        unifac_result={"DOH": 1},
        psrk_result={"DOH": 1},
        joback_result={"-CH2-": 2, "-OH (alcohol)": 2},
    ),
    Case(
        identifier="CCOCCO",
        identifier_type="smiles",
        cases_module="alcohols",
        r=None,
        q=None,
        unifac_result={"CH3": 1, "CH2": 1, "C2H5O2": 1},
        psrk_result={"CH3": 1, "CH2": 1, "C2H5O2": 1},
        joback_result={
            "-CH3": 1,
            "-CH2-": 3,
            "-OH (alcohol)": 1,
            "-O- (non-ring)": 1,
        },
    ),
    Case(
        identifier="CCOC(C)CO",
        identifier_type="smiles",
        cases_module="alcohols",
        r=None,
        q=None,
        unifac_result={"CH3": 2, "C2H4O2": 1, "CH2": 1},
        psrk_result={"CH3": 2, "C2H4O2": 1, "CH2": 1},
        joback_result={
            "-CH3": 2,
            "-CH2-": 2,
            ">CH-": 1,
            "-OH (alcohol)": 1,
            "-O- (non-ring)": 1,
        },
    ),
    Case(
        "CC(C)(C)C1CCC(C(C1)O)(C2=CC=CC=C2)O",
        "smiles",
        "alcohols",
        "1,2-Cyclohexanediol, 4-tert-butyl-1-phenyl-, stereoisomer",
        r=None,
        q=None,
        unifac_result={
            "CH3": 3,
            "CH2": 3,
            "CH": 2,
            "C": 2,
            "OH": 2,
            "AC": 1,
            "ACH": 5,
        },
        psrk_result={
            "CH3": 3,
            "CH2": 3,
            "CH": 2,
            "C": 2,
            "OH": 2,
            "AC": 1,
            "ACH": 5,
        },
        joback_result={
            "-CH3": 3,
            ">C<": 1,
            "ring-CH2-": 3,
            "ring>CH-": 2,
            "ring>C<": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-OH (alcohol)": 2,
        },
    ),
    Case(
        "CCCC(C(C)CO)O",
        "smiles",
        "alcohols",
        "(2S,3S)-2-Methyl-1,3-hexanediol",
        r=None,
        q=None,
        unifac_result={"CH3": 2, "CH2": 3, "CH": 2, "OH": 2},
        psrk_result={"CH3": 2, "CH2": 3, "CH": 2, "OH": 2},
        joback_result={"-CH3": 2, "-CH2-": 3, ">CH-": 2, "-OH (alcohol)": 2},
    ),
    Case(
        "CC(C)O",
        "smiles",
        "alcohols",
        "2-propanol",
        r=None,
        q=None,
        unifac_result={"CH3": 2, "CH": 1, "OH": 1},
        psrk_result={"CH3": 2, "CH": 1, "OH": 1},
        joback_result={"-CH3": 2, ">CH-": 1, "-OH (alcohol)": 1},
    ),
    Case(
        "CO",
        "smiles",
        "alcohols",
        "methanol",
        r=None,
        q=None,
        unifac_result={"CH3OH": 1},
        psrk_result={"CH3OH": 1},
        joback_result={"-CH3": 1, "-OH (alcohol)": 1},
    ),
    Case(
        "CCO",
        "smiles",
        "alcohols",
        "ethanol",
        r=None,
        q=None,
        unifac_result={"CH3": 1, "CH2": 1, "OH": 1},
        psrk_result={"CH3": 1, "CH2": 1, "OH": 1},
        joback_result={"-CH3": 1, "-CH2-": 1, "-OH (alcohol)": 1},
    ),
    Case(
        "CCCO",
        "smiles",
        "alcohols",
        "1-propanol",
        r=None,
        q=None,
        unifac_result={"CH3": 1, "CH2": 2, "OH": 1},
        psrk_result={"CH3": 1, "CH2": 2, "OH": 1},
        joback_result={"-CH3": 1, "-CH2-": 2, "-OH (alcohol)": 1},
    ),
    Case(
        "C1=CC=C2C(=C1)C=CC3=C2C(=C(C=C3)O)O",
        "smiles",
        "alcohols",
        "Phenanthrene-3,4-diol",
        r=None,
        q=None,
        unifac_result={"ACH": 8, "AC": 4, "ACOH": 2},
        psrk_result={"ACH": 8, "AC": 4, "ACOH": 2},
        joback_result={"ring=CH-": 8, "ring=C<": 6, "-OH (phenol)": 2},
    ),
    Case(
        "CC(C)(C)C1=C(C(=CC=C1)O)O",
        "smiles",
        "alcohols",
        "3-(tert-butyl)benzene-1,2-diol",
        r=None,
        q=None,
        unifac_result={"ACH": 3, "AC": 1, "ACOH": 2, "CH3": 3, "C": 1},
        psrk_result={"ACH": 3, "AC": 1, "ACOH": 2, "CH3": 3, "C": 1},
        joback_result={
            "-CH3": 3,
            ">C<": 1,
            "ring=CH-": 3,
            "ring=C<": 3,
            "-OH (phenol)": 2,
        },
    ),
    Case(
        "C1=CC(=CC(=C1)O)C2=C(C=C(C=C2)O)O",
        "smiles",
        "alcohols",
        "[1,1'-Biphenyl]-2,3',4-triol",
        r=None,
        q=None,
        unifac_result={"ACH": 7, "AC": 2, "ACOH": 3},
        psrk_result={"ACH": 7, "AC": 2, "ACOH": 3},
        joback_result={"ring=CH-": 7, "ring=C<": 5, "-OH (phenol)": 3},
    ),
]
