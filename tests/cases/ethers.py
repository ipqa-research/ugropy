# =============================================================================
# Description: This file contains the test cases for the ethers molecules.
#
# There is a lot of combinations for ethers molecules, also are the most
# problematics. TODO: moar tests
# =============================================================================
from .case import Case


ethers_cases = [
    Case(
        identifier="OC1CC(OC2=CC=CC=C12)C1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={
            "ACH": 9,
            "AC": 2,
            "ACCH": 1,
            "CHO": 1,
            "CH2": 1,
            "OH": 1,
        },
        psrk_result={
            "ACH": 9,
            "AC": 2,
            "ACCH": 1,
            "CHO": 1,
            "CH2": 1,
            "OH": 1,
        },
        joback_result={
            "ring-CH2-": 1,
            "ring>CH-": 2,
            "ring=CH-": 9,
            "ring=C<": 3,
            "-OH (alcohol)": 1,
            "-O- (ring)": 1,
        },
        dortmund_result={
            "ACH": 9,
            "AC": 3,
            "CY-CH": 1,
            "CHO": 1,
            "CY-CH2": 1,
            "OH (S)": 1,
        },
    ),
    Case(
        identifier="O[C@@H]1CO[C@H](O)[C@@H](O)[C@@H]1O",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH2": 1, "CH": 3, "OH": 4, "CHO": 1},
            {"CH": 4, "OH": 4, "CH2O": 1},
        ],
        psrk_result=[
            {"CH2": 1, "CH": 3, "OH": 4, "CHO": 1},
            {"CH": 4, "OH": 4, "CH2O": 1},
        ],
        joback_result={
            "ring-CH2-": 1,
            "ring>CH-": 4,
            "-OH (alcohol)": 4,
            "-O- (ring)": 1,
        },
        dortmund_result=[
            {"OH (P)": 1, "CHO": 1, "CY-CH2": 1, "CY-CH": 3, "OH (S)": 3},
            {"OH (P)": 1, "CY-CH": 4, "OH (S)": 3, "CY-CH2O": 1},
        ],
    ),
    Case(
        identifier="C1COCCOCCOCCOC1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH2O": 4, "CH2": 5},
        psrk_result={"CH2O": 4, "CH2": 5},
        joback_result={"ring-CH2-": 9, "-O- (ring)": 4},
        dortmund_result={"CY-CH2O": 4, "CY-CH2": 5},
    ),
    Case(
        identifier="C1COCCO1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH2O": 2, "CH2": 2},
        psrk_result={"CH2O": 2, "CH2": 2},
        joback_result={"ring-CH2-": 4, "-O- (ring)": 2},
        dortmund_result={"CY-CH2O": 2, "CY-CH2": 2},
    ),
    Case(
        identifier="C1COCO1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH2O": 2, "CH2": 1},
        psrk_result={"CH2O": 2, "CH2": 1},
        joback_result={"ring-CH2-": 3, "-O- (ring)": 2},
        dortmund_result={"CY-CH2O": 2, "CY-CH2": 1},
    ),
    Case(
        identifier="C1COCCOCCOCCOCCOCCO1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH2O": 6, "CH2": 6},
        psrk_result={"CH2O": 6, "CH2": 6},
        joback_result={"ring-CH2-": 12, "-O- (ring)": 6},
        dortmund_result={"CY-CH2O": 6, "CY-CH2": 6},
    ),
    Case(
        identifier="C1CCOC1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"THF": 1, "CH2": 3},
        psrk_result={"THF": 1, "CH2": 3},
        joback_result={"ring-CH2-": 4, "-O- (ring)": 1},
        dortmund_result={"THF": 1, "CY-CH2": 2},
    ),
    Case(
        identifier="CC1COCC1C",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"THF": 1, "CH2": 1, "CH": 2, "CH3": 2},
        psrk_result={"THF": 1, "CH2": 1, "CH": 2, "CH3": 2},
        joback_result={
            "-CH3": 2,
            "ring-CH2-": 2,
            "ring>CH-": 2,
            "-O- (ring)": 1,
        },
        dortmund_result={"THF": 1, "CY-CH": 2, "CH3": 2},
    ),
    Case(
        identifier="CC1COCC1O",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"THF": 1, "CH2": 1, "CH": 2, "CH3": 1, "OH": 1},
        psrk_result={"THF": 1, "CH2": 1, "CH": 2, "CH3": 1, "OH": 1},
        joback_result={
            "-CH3": 1,
            "ring-CH2-": 2,
            "ring>CH-": 2,
            "-OH (alcohol)": 1,
            "-O- (ring)": 1,
        },
        dortmund_result={"THF": 1, "CY-CH": 2, "CH3": 1, "OH (S)": 1},
    ),
    Case(
        identifier="CC(C)OC(C)C",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 4, "CH": 1, "CHO": 1},
        psrk_result={"CH3": 4, "CH": 1, "CHO": 1},
        joback_result={"-CH3": 4, ">CH-": 2, "-O- (non-ring)": 1},
        dortmund_result={"CH3": 4, "CH": 1, "CHO": 1},
    ),
    Case(
        identifier="CCOCC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 2, "CH2": 1, "CH2O": 1},
        psrk_result={"CH3": 2, "CH2": 1, "CH2O": 1},
        joback_result={"-CH3": 2, "-CH2-": 2, "-O- (non-ring)": 1},
        dortmund_result={"CH3": 2, "CH2": 1, "CH2O": 1},
    ),
    Case(
        identifier="COC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 1, "CH3O": 1},
        psrk_result={"CH3": 1, "CH3O": 1},
        joback_result={"-CH3": 2, "-O- (non-ring)": 1},
        dortmund_result={"CH3": 1, "CH3O": 1},
    ),
    Case(
        identifier="C1CCC(CC1)OC2CCCCO2",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH2": 8, "CH": 1, "CH2O": 1, "CHO": 1},
            {"CH2": 9, "CHO": 2},
        ],
        psrk_result=[
            {"CH2": 8, "CH": 1, "CH2O": 1, "CHO": 1},
            {"CH2": 9, "CHO": 2},
        ],
        joback_result={
            "ring-CH2-": 9,
            "ring>CH-": 2,
            "-O- (non-ring)": 1,
            "-O- (ring)": 1,
        },
        dortmund_result=[
            {"CHO": 1, "CY-CH2": 8, "CY-CH": 1, "CY-CH2O": 1},
            {"CHO": 2, "CY-CH2": 9},
        ],
    ),
    Case(
        identifier="CCOC(=O)OC1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 1, "CH2O": 1, "COO": 1, "AC": 1, "ACH": 5},
        psrk_result={"CH3": 1, "CH2O": 1, "COO": 1, "AC": 1, "ACH": 5},
        joback_result={
            "-CH3": 1,
            "-CH2-": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result={"CH3": 1, "CH2O": 1, "COO": 1, "AC": 1, "ACH": 5},
    ),
    Case(
        identifier="CC(C)OC(=O)OC1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 2, "CHO": 1, "COO": 1, "AC": 1, "ACH": 5},
        psrk_result={"CH3": 2, "CHO": 1, "COO": 1, "AC": 1, "ACH": 5},
        joback_result={
            "-CH3": 2,
            ">CH-": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result={"CH3": 2, "CHO": 1, "COO": 1, "AC": 1, "ACH": 5},
    ),
    Case(
        identifier="CC(C)(C)OC(=O)OC1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={},
        psrk_result={},
        joback_result={
            "-CH3": 3,
            ">C<": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result={},
    ),
    Case(
        identifier="C1=CC=C(C=C1)COC(=O)OCCO",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"ACCH2": 1, "ACH": 5, "COO": 1, "C2H5O2": 1},
        psrk_result={"ACCH2": 1, "ACH": 5, "COO": 1, "C2H5O2": 1},
        joback_result={
            "-CH2-": 3,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-OH (alcohol)": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result={
            "CH2": 1,
            "ACH": 5,
            "AC": 1,
            "OH (P)": 1,
            "(CH2)2CB": 1,
        },
    ),
    Case(
        identifier="CCOC(=O)OC(C)(C)C",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 4, "C": 1, "COO": 1, "CH2O": 1},
        psrk_result={"CH3": 4, "C": 1, "COO": 1, "CH2O": 1},
        joback_result={
            "-CH3": 4,
            "-CH2-": 1,
            ">C<": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result={"CH3": 4, "C": 1, "COO": 1, "CH2O": 1},
    ),
    Case(
        identifier="CCOC(=O)OC1C(CCC(C1C)C)C",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH3": 4, "CH2": 2, "CH": 4, "CH2O": 1, "COO": 1},
            {"CH3": 4, "CH2": 3, "CH": 3, "CHO": 1, "COO": 1},
        ],
        psrk_result=[
            {"CH3": 4, "CH2": 2, "CH": 4, "CH2O": 1, "COO": 1},
            {"CH3": 4, "CH2": 3, "CH": 3, "CHO": 1, "COO": 1},
        ],
        joback_result={
            "-CH3": 4,
            "-CH2-": 1,
            "ring-CH2-": 2,
            "ring>CH-": 4,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result=[
            {"CH3": 4, "CH2": 1, "CHO": 1, "COO": 1, "CY-CH2": 2, "CY-CH": 3},
            {"CH3": 4, "CH2O": 1, "COO": 1, "CY-CH2": 2, "CY-CH": 4},
        ],
    ),
    Case(
        identifier="CCOC(=O)OCC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 2, "CH2": 1, "COO": 1, "CH2O": 1},
        psrk_result={"CH3": 2, "CH2": 1, "COO": 1, "CH2O": 1},
        joback_result={
            "-CH3": 2,
            "-CH2-": 2,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result={"CH3": 2, "(CH2)2CB": 1},
    ),
    Case(
        identifier="COC(=O)OC1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"AC": 1, "ACH": 5, "COO": 1, "CH3O": 1},
        psrk_result={"AC": 1, "ACH": 5, "COO": 1, "CH3O": 1},
        joback_result={
            "-CH3": 1,
            "ring=CH-": 5,
            "ring=C<": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result={"AC": 1, "ACH": 5, "COO": 1, "CH3O": 1},
    ),
    Case(
        identifier="CC(C)(C)OC(=O)OC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 3, "C": 1, "COO": 1, "CH3O": 1},
        psrk_result={"CH3": 3, "C": 1, "COO": 1, "CH3O": 1},
        joback_result={
            "-CH3": 4,
            ">C<": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result={"CH3": 3, "C": 1, "COO": 1, "CH3O": 1},
    ),
    Case(
        identifier="CC(C)OC(=O)OC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH3": 3, "CHO": 1, "COO": 1},
            {"CH3": 2, "CH": 1, "CH3O": 1, "COO": 1},
        ],
        psrk_result=[
            {"CH3": 3, "CHO": 1, "COO": 1},
            {"CH3": 2, "CH": 1, "CH3O": 1, "COO": 1},
        ],
        joback_result={
            "-CH3": 3,
            ">CH-": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result=[
            {"CH3": 3, "CHO": 1, "COO": 1},
            {"CH3": 2, "CH": 1, "CH3O": 1, "COO": 1},
        ],
    ),
    Case(
        identifier="CCOC(=O)OC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH3": 2, "CH2O": 1, "COO": 1},
            {"CH3": 1, "CH2": 1, "CH3O": 1, "COO": 1},
        ],
        psrk_result=[
            {"CH3": 2, "CH2O": 1, "COO": 1},
            {"CH3": 1, "CH2": 1, "CH3O": 1, "COO": 1},
        ],
        joback_result={
            "-CH3": 2,
            "-CH2-": 1,
            "-O- (non-ring)": 1,
            "-COO- (ester)": 1,
        },
        dortmund_result={"CH3": 1, "CH2CH3CB": 1},
    ),
    Case(
        identifier="COC(=O)OC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 1, "COO": 1, "CH3O": 1},
        psrk_result={"CH3": 1, "COO": 1, "CH3O": 1},
        joback_result={"-CH3": 2, "-O- (non-ring)": 1, "-COO- (ester)": 1},
        dortmund_result={"(CH3)2CB": 1},
    ),
    Case(
        identifier="COCOC(C)OCOC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH3": 1, "CH": 1, "CH3O": 2, "CH2O": 2},
            {"CH3": 2, "CH3O": 1, "CH2O": 2, "CHO": 1},
            {"CH3": 1, "CH2": 1, "CH3O": 2, "CH2O": 1, "CHO": 1},
        ],
        psrk_result=[
            {"CH3": 1, "CH": 1, "CH3O": 2, "CH2O": 2},
            {"CH3": 2, "CH3O": 1, "CH2O": 2, "CHO": 1},
            {"CH3": 1, "CH2": 1, "CH3O": 2, "CH2O": 1, "CHO": 1},
        ],
        joback_result={"-CH3": 3, "-CH2-": 2, ">CH-": 1, "-O- (non-ring)": 4},
        dortmund_result=[
            {"CH3": 1, "CH": 1, "CH3O": 2, "CH2O": 2},
            {"CH3": 2, "CH3O": 1, "CH2O": 2, "CHO": 1},
            {"CH3": 1, "CH2": 1, "CH3O": 2, "CH2O": 1, "CHO": 1},
        ],
    ),
    Case(
        identifier="CC(C)OCOC(C)OCOC(C)C",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH3": 5, "CH": 1, "CHO": 2, "CH2O": 2},
            {"CH3": 5, "CH2": 1, "CH2O": 1, "CHO": 3},
        ],
        psrk_result=[
            {"CH3": 5, "CH": 1, "CHO": 2, "CH2O": 2},
            {"CH3": 5, "CH2": 1, "CH2O": 1, "CHO": 3},
        ],
        joback_result={"-CH3": 5, "-CH2-": 2, ">CH-": 3, "-O- (non-ring)": 4},
        dortmund_result=[
            {"CH3": 5, "CH": 1, "CHO": 2, "CH2O": 2},
            {"CH3": 5, "CH2": 1, "CH2O": 1, "CHO": 3},
        ],
    ),
    Case(
        identifier="CC(C)OCOCC(OCOC(C)C)OCOC(C)C",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH3": 6, "CH": 2, "CH2O": 4, "CHO": 2},
            {"CH3": 6, "CH2": 1, "CH": 1, "CH2O": 3, "CHO": 3},
            {"CH3": 6, "CH2": 2, "CH2O": 2, "CHO": 4},
        ],
        psrk_result=[
            {"CH3": 6, "CH": 2, "CH2O": 4, "CHO": 2},
            {"CH3": 6, "CH2": 1, "CH": 1, "CH2O": 3, "CHO": 3},
            {"CH3": 6, "CH2": 2, "CH2O": 2, "CHO": 4},
        ],
        joback_result={"-CH3": 6, "-CH2-": 4, ">CH-": 4, "-O- (non-ring)": 6},
        dortmund_result=[
            {"CH3": 6, "CH": 2, "CH2O": 4, "CHO": 2},
            {"CH3": 6, "CH2": 1, "CH": 1, "CH2O": 3, "CHO": 3},
            {"CH3": 6, "CH2": 2, "CH2O": 2, "CHO": 4},
        ],
    ),
    Case(
        identifier="CC(C)OCOC(OCOC(C)C)OCOC(C)C",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH3": 6, "CHO": 3, "CH2O": 3, "CH": 1},
            {"CH3": 6, "CH2": 1, "CH2O": 2, "CHO": 4},
        ],
        psrk_result=[
            {"CH3": 6, "CHO": 3, "CH2O": 3, "CH": 1},
            {"CH3": 6, "CH2": 1, "CH2O": 2, "CHO": 4},
        ],
        joback_result={"-CH3": 6, "-CH2-": 3, ">CH-": 4, "-O- (non-ring)": 6},
        dortmund_result=[
            {"CH3": 6, "CHO": 3, "CH2O": 3, "CH": 1},
            {"CH3": 6, "CH2": 1, "CH2O": 2, "CHO": 4},
        ],
    ),
    Case(
        identifier="CC(C)OCOC(C)C",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH3": 4, "CHO": 1, "CH2O": 1, "CH": 1},
            {"CH3": 4, "CH2": 1, "CHO": 2},
        ],
        psrk_result=[
            {"CH3": 4, "CHO": 1, "CH2O": 1, "CH": 1},
            {"CH3": 4, "CH2": 1, "CHO": 2},
        ],
        joback_result={"-CH3": 4, "-CH2-": 1, ">CH-": 2, "-O- (non-ring)": 2},
        dortmund_result=[
            {"CH3": 4, "CHO": 1, "CH2O": 1, "CH": 1},
            {"CH3": 4, "CH2": 1, "CHO": 2},
        ],
    ),
    Case(
        identifier="CCOCOCC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH3": 2, "CH2O": 2, "CH2": 1},
        psrk_result={"CH3": 2, "CH2O": 2, "CH2": 1},
        joback_result={"-CH3": 2, "-CH2-": 3, "-O- (non-ring)": 2},
        dortmund_result={"CH3": 2, "CH2O": 2, "CH2": 1},
    ),
    Case(
        identifier="COCOC",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result=[
            {"CH3O": 2, "CH2": 1},
            {"CH3": 1, "CH3O": 1, "CH2O": 1},
        ],
        psrk_result=[
            {"CH3O": 2, "CH2": 1},
            {"CH3": 1, "CH3O": 1, "CH2O": 1},
        ],
        joback_result={"-CH3": 2, "-CH2-": 1, "-O- (non-ring)": 2},
        dortmund_result=[
            {"CH3O": 2, "CH2": 1},
            {"CH3": 1, "CH3O": 1, "CH2O": 1},
        ],
    ),
    Case(
        identifier="COC(O)=O",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"COOH": 1, "CH3O": 1},
        psrk_result={"COOH": 1, "CH3O": 1},
        joback_result=[
            {"-CH3": 1, "-O- (non-ring)": 1, "-COOH (acid)": 1},
            {"-CH3": 1, "-OH (alcohol)": 1, "-COO- (ester)": 1},
        ],
        dortmund_result={"COOH": 1, "CH3O": 1},
    ),
    Case(
        identifier="CC(C)OC(O)=O",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"COOH": 1, "CHO": 1, "CH3": 2},
        psrk_result={"COOH": 1, "CHO": 1, "CH3": 2},
        joback_result=[
            {"-CH3": 2, ">CH-": 1, "-O- (non-ring)": 1, "-COOH (acid)": 1},
            {"-CH3": 2, ">CH-": 1, "-OH (alcohol)": 1, "-COO- (ester)": 1},
        ],
        dortmund_result={"COOH": 1, "CHO": 1, "CH3": 2},
    ),
    Case(
        identifier="CC(C)(C)OC(O)=O",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"OH": 1, "COO": 1, "C": 1, "CH3": 3},
        psrk_result={"OH": 1, "COO": 1, "C": 1, "CH3": 3},
        joback_result=[
            {"-CH3": 3, ">C<": 1, "-O- (non-ring)": 1, "-COOH (acid)": 1},
            {"-CH3": 3, ">C<": 1, "-OH (alcohol)": 1, "-COO- (ester)": 1},
        ],
        dortmund_result={},
    ),
    Case(
        identifier="OC(=O)OC1=CC=CC=C1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"OH": 1, "COO": 1, "AC": 1, "ACH": 5},
        psrk_result={"OH": 1, "COO": 1, "AC": 1, "ACH": 5},
        joback_result=[
            {
                "ring=CH-": 5,
                "ring=C<": 1,
                "-OH (alcohol)": 1,
                "-COO- (ester)": 1,
            },
            {
                "ring=CH-": 5,
                "ring=C<": 1,
                "-O- (non-ring)": 1,
                "-COOH (acid)": 1,
            },
        ],
        dortmund_result={},
    ),
    Case(
        identifier="C1COCON1",
        identifier_type="smiles",
        cases_module="ethers",
        unifac_result={"CH2O": 2, "CH2NH": 1},
        psrk_result={"CH2O": 2, "CH2NH": 1},
        joback_result={"ring-CH2-": 3, "-O- (ring)": 2, ">NH (ring)": 1},
        dortmund_result={"CY-CH2O": 2, "CH2NH": 1},
    ),
]
