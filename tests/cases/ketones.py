# =============================================================================
# Description: This file contains the test cases for the ketones molecules
#
# If possible simple ketones, don't think too much about it, if your objetive
# is a ketone detection put it here.
# =============================================================================
from .case import Case


ketones_cases = [
    Case(
        "CC(C)CC(=O)OCC(=O)C(C)O",
        "smiles",
        "ketones",
        "2,2-Dimethyl-1,3-propanediol diacetate",
        unifac_result={"CH3": 3, "CH": 2, "CH2COO": 1, "CH2CO": 1, "OH": 1},
        psrk_result={"CH3": 3, "CH": 2, "CH2COO": 1, "CH2CO": 1, "OH": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 2,
            ">CH-": 2,
            "-OH (alcohol)": 1,
            ">C=O (non-ring)": 1,
            "-COO- (ester)": 1,
        },
    ),
    Case(
        "O=C(CC1=CC=CC=C1)CC1=CC=CC=C1",
        "smiles",
        "ketones",
        unifac_result={"ACH": 10, "AC": 1, "ACCH2": 1, "CH2CO": 1},
        psrk_result={"ACH": 10, "AC": 1, "ACCH2": 1, "CH2CO": 1},
        joback_result={
            "-CH2-": 2,
            "ring=CH-": 10,
            "ring=C<": 2,
            ">C=O (non-ring)": 1,
        },
    ),
    Case(
        "CC(=O)CC1=CC=CC=C1",
        "smiles",
        "ketones",
        unifac_result={"CH3CO": 1, "ACH": 5, "ACCH2": 1},
        psrk_result={"CH3CO": 1, "ACH": 5, "ACCH2": 1},
        joback_result={
            "-CH3": 1,
            "-CH2-": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            ">C=O (non-ring)": 1,
        },
    ),
    Case(
        "CC(C)(C)C(=O)CC1=CC=CC=C1",
        "smiles",
        "ketones",
        unifac_result={"CH3": 3, "ACH": 5, "AC": 1, "CH2CO": 1, "C": 1},
        psrk_result={"CH3": 3, "ACH": 5, "AC": 1, "CH2CO": 1, "C": 1},
        joback_result={
            "-CH3": 3,
            "-CH2-": 1,
            ">C<": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            ">C=O (non-ring)": 1,
        },
    ),
    Case(
        "C1CC1=O",
        "smiles",
        "ketones",
        commentary="Cyclopropanone",
        unifac_result={"CH2": 1, "CH2CO": 1},
        psrk_result={"CH2": 1, "CH2CO": 1},
        joback_result={"ring-CH2-": 2, ">C=O (ring)": 1},
    ),
    Case(
        "C1CCCC=CCCCCCCCC(=O)CCC1",
        "smiles",
        "ketones",
        commentary="(9Z)-Cycloheptadec-9-en-1-one",
        unifac_result={"CH2": 13, "CH=CH": 1, "CH2CO": 1},
        psrk_result={"CH2": 13, "CH=CH": 1, "CH2CO": 1},
        joback_result={"ring-CH2-": 14, "ring=CH-": 2, ">C=O (ring)": 1},
    ),
    Case(
        "CC1=CC(=CC(=C1C(=O)C)C)C(C)(C)C",
        "smiles",
        "ketones",
        commentary="1-(4-tert-butyl-2,6-dimethylphenyl)ethanone",
        unifac_result={
            "ACH": 2,
            "AC": 2,
            "ACCH3": 2,
            "CH3": 3,
            "C": 1,
            "CH3CO": 1,
        },
        psrk_result={
            "ACH": 2,
            "AC": 2,
            "ACCH3": 2,
            "CH3": 3,
            "C": 1,
            "CH3CO": 1,
        },
        joback_result={
            "-CH3": 6,
            ">C<": 1,
            "ring=CH-": 2,
            "ring=C<": 4,
            ">C=O (non-ring)": 1,
        },
    ),
    Case(
        "CC(=O)C1=CC=CC=C1",
        "smiles",
        "ketones",
        commentary="Acetophenone",
        unifac_result={"ACH": 5, "AC": 1, "CH3CO": 1},
        psrk_result={"ACH": 5, "AC": 1, "CH3CO": 1},
        joback_result={
            "-CH3": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            ">C=O (non-ring)": 1,
        },
    ),
    Case(
        "CC(=O)C",
        "smiles",
        "ketones",
        commentary="acetone",
        unifac_result={"CH3CO": 1, "CH3": 1},
        psrk_result={"CH3CO": 1, "CH3": 1},
        joback_result={"-CH3": 2, ">C=O (non-ring)": 1},
    ),
    Case(
        "CCC(=O)CC",
        "smiles",
        "ketones",
        commentary="3-pentanone",
        unifac_result={"CH3": 2, "CH2": 1, "CH2CO": 1},
        psrk_result={"CH3": 2, "CH2": 1, "CH2CO": 1},
        joback_result={"-CH3": 2, "-CH2-": 2, ">C=O (non-ring)": 1},
    ),
    Case(
        "CCC(=O)C",
        "smiles",
        "ketones",
        commentary="2-butanone",
        unifac_result=[
            {"CH3": 1, "CH2": 1, "CH3CO": 1},
            {"CH3": 2, "CH2CO": 1},
        ],
        psrk_result=[
            {"CH3": 1, "CH2": 1, "CH3CO": 1},
            {"CH3": 2, "CH2CO": 1},
        ],
        joback_result={"-CH3": 2, "-CH2-": 1, ">C=O (non-ring)": 1},
    ),
]
