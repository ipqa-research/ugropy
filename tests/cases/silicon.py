# =============================================================================
# Silicon test cases
# =============================================================================
from .case import Case


silicon_cases = [
    Case(
        identifier="C[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 1, "SIH3": 1},
        psrk_result={"CH3": 1, "SIH3": 1},
        joback_result={},
    ),
    Case(
        identifier="CC[Si](CC)([H])[H]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 2, "CH2": 2, "SIH2": 1},
        psrk_result={"CH3": 2, "CH2": 2, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="C[Si](O[Si](C)(C)C)(O[Si](C)(C)C)[H]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 7, "SI": 1, "SIHO": 1, "SIO": 1},
            {"CH3": 7, "SIH": 1, "SIO": 2},
        ],
        psrk_result={"CH3": 7, "SI": 1, "SIHO": 1, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="C[Si](C)(C)O[Si](C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 6, "SI": 1, "SIO": 1},
        psrk_result={"CH3": 6, "SI": 1, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH2]O[SiH](C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "SIH2": 1, "SIHO": 1},
            {"CH3": 3, "SIH": 1, "SIH2O": 1},
        ],
        psrk_result={"CH3": 3, "SIH": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH2]O[Si](C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 4, "SIH2": 1, "SIO": 1},
            {"CH3": 4, "SI": 1, "SIH2O": 1},
        ],
        psrk_result={"CH3": 4, "SI": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH](C)O[Si](C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 5, "SIH": 1, "SIO": 1},
            {"CH3": 5, "SI": 1, "SIHO": 1},
        ],
        psrk_result={"CH3": 5, "SI": 1, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)(C)O[Si](C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 6, "C": 1, "SIO": 1},
        psrk_result={"CH3": 6, "C": 1, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)O[Si](C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 5, "CH": 1, "SIO": 1},
            {"CH3": 5, "CHO": 1, "SI": 1},
        ],
        psrk_result={"CH3": 5, "CHO": 1, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)O[SiH](C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 4, "CH": 1, "SIHO": 1},
            {"CH3": 4, "CHO": 1, "SIH": 1},
        ],
        psrk_result={"CH3": 4, "CHO": 1, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH2]OC(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CH": 1, "SIH2O": 1},
            {"CH3": 3, "CHO": 1, "SIH2": 1},
        ],
        psrk_result={"CH3": 3, "CHO": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="CCO[SiH2]C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "CH2": 1, "SIH2O": 1},
            {"CH3": 2, "CH2O": 1, "SIH2": 1},
        ],
        psrk_result={"CH3": 2, "CH2O": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="CO[SiH2]C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "SIH2O": 1},
            {"CH3": 1, "CH3O": 1, "SIH2": 1},
        ],
        psrk_result={"CH3": 1, "CH3O": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="C[Si](O[Si](C)([H])[H])([H])[H]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 2, "SIH2": 1, "SIH2O": 1},
        psrk_result={"CH3": 2, "SIH2": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="C[Si](C)(O[Si](C)(C)[H])[H]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 4, "SIH": 1, "SIHO": 1},
        psrk_result={"CH3": 4, "SIH": 1, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="C[Si]1(O[Si](O[Si](O[Si](O1)(C)C)(C)C)(C)C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 8, "SIO": 4},
        psrk_result={"CH3": 8, "SIO": 4},
        joback_result={},
    ),
    Case(
        identifier="CC(=O)O[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3COO": 1, "SIH3": 1},
        psrk_result={"CH3COO": 1, "SIH3": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(=O)O[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3COO": 1, "SIH3": 1, "SIH2": 1},
        psrk_result={"CH3COO": 1, "SIH3": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(=O)O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3COO": 1, "SIH3": 2, "SIH": 1},
        psrk_result={"CH3COO": 1, "SIH3": 2, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3COO": 1, "SIH3": 3, "SI": 1},
        psrk_result={"CH3COO": 1, "SIH3": 3, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="CCC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 1, "CH2COO": 1, "SIH3": 3, "SI": 1},
        psrk_result={"CH3": 1, "CH2COO": 1, "SIH3": 3, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)C(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 2, "CH": 1, "COO": 1, "SIH3": 3, "SI": 1},
        psrk_result={"CH3": 2, "CH": 1, "COO": 1, "SIH3": 3, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3]OC(=O)O[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"COO": 1, "SIH3": 2, "SIH2O": 1},
        psrk_result={"COO": 1, "SIH3": 2, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH2]OC(=O)O[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"COO": 1, "SIH3": 2, "SIH2O": 1, "SIH2": 1},
        psrk_result={"COO": 1, "SIH3": 2, "SIH2O": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH2]OC(=O)O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"COO": 1, "SIH3": 3, "SIH2": 1, "SIHO": 1},
            {"COO": 1, "SIH3": 3, "SIH": 1, "SIH2O": 1},
        ],
        psrk_result={"COO": 1, "SIH3": 3, "SIH2O": 1, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH2]OC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"COO": 1, "SIH3": 4, "SIH2": 1, "SIO": 1},
            {"COO": 1, "SIH3": 4, "SI": 1, "SIH2O": 1},
        ],
        psrk_result={"COO": 1, "SIH3": 4, "SIH2O": 1, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3]OC(=O)O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"COO": 1, "SIH3": 3, "SIHO": 1},
        psrk_result={"COO": 1, "SIH3": 3, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH]([SiH3])OC(=O)O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"COO": 1, "SIH3": 4, "SIHO": 1, "SIH": 1},
        psrk_result={"COO": 1, "SIH3": 4, "SIHO": 1, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH]([SiH3])OC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"COO": 1, "SIH3": 5, "SIH": 1, "SIO": 1},
            {"COO": 1, "SIH3": 5, "SI": 1, "SIHO": 1},
        ],
        psrk_result={"COO": 1, "SIH3": 5, "SIHO": 1, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3]OC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"COO": 1, "SIH3": 4, "SIO": 1},
        psrk_result={"COO": 1, "SIH3": 4, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][Si]([SiH3])([SiH3])OC(=O)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"COO": 1, "SIH3": 6, "SIO": 1, "SI": 1},
        psrk_result={"COO": 1, "SIH3": 6, "SIO": 1, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3]OC=O",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"HCOO": 1, "SIH3": 1},
        psrk_result={"HCOO": 1, "SIH3": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH2]OC=O",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"HCOO": 1, "SIH3": 1, "SIH2": 1},
        psrk_result={"HCOO": 1, "SIH3": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH]([SiH3])OC=O",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"HCOO": 1, "SIH3": 2, "SIH": 1},
        psrk_result={"HCOO": 1, "SIH3": 2, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][Si]([SiH3])([SiH3])OC=O",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"HCOO": 1, "SIH3": 3, "SI": 1},
        psrk_result={"HCOO": 1, "SIH3": 3, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="COC(=O)O[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3O": 1, "COO": 1, "SIH3": 1},
        psrk_result={"CH3O": 1, "COO": 1, "SIH3": 1},
        joback_result={},
    ),
    Case(
        identifier="COC(=O)O[SiH2]C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "COO": 1, "SIH2O": 1},
            {"CH3": 1, "CH3O": 1, "COO": 1, "SIH2": 1},
        ],
        psrk_result={"CH3": 1, "CH3O": 1, "COO": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="COC(=O)O[SiH](C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "COO": 1, "SIHO": 1},
            {"CH3": 2, "CH3O": 1, "COO": 1, "SIH": 1},
        ],
        psrk_result={"CH3": 2, "CH3O": 1, "COO": 1, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="COC(=O)O[Si](C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 4, "COO": 1, "SIO": 1},
            {"CH3": 3, "CH3O": 1, "COO": 1, "SI": 1},
        ],
        psrk_result={"CH3": 3, "CH3O": 1, "COO": 1, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="CCOC(=O)O[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 1, "CH2O": 1, "COO": 1, "SIH3": 1},
        psrk_result={"CH3": 1, "CH2O": 1, "COO": 1, "SIH3": 1},
        joback_result={},
    ),
    Case(
        identifier="CCOC(=O)O[SiH2]C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "CH2": 1, "COO": 1, "SIH2O": 1},
            {"CH3": 2, "CH2O": 1, "COO": 1, "SIH2": 1},
        ],
        psrk_result={"CH3": 2, "CH2O": 1, "COO": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="CCOC(=O)O[SiH](C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CH2": 1, "COO": 1, "SIHO": 1},
            {"CH3": 3, "CH2O": 1, "COO": 1, "SIH": 1},
        ],
        psrk_result={"CH3": 3, "CH2O": 1, "COO": 1, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="CCOC(=O)O[Si](C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 4, "CH2": 1, "COO": 1, "SIO": 1},
            {"CH3": 4, "CH2O": 1, "COO": 1, "SI": 1},
        ],
        psrk_result={"CH3": 4, "CH2O": 1, "COO": 1, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)OC(=O)O[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 2, "CHO": 1, "COO": 1, "SIH3": 1},
        psrk_result={"CH3": 2, "CHO": 1, "COO": 1, "SIH3": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH2]OC(=O)OC(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CH": 1, "COO": 1, "SIH2O": 1},
            {"CH3": 3, "CHO": 1, "COO": 1, "SIH2": 1},
        ],
        psrk_result={"CH3": 3, "CHO": 1, "COO": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)OC(=O)O[SiH](C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 4, "CH": 1, "COO": 1, "SIHO": 1},
            {"CH3": 4, "CHO": 1, "COO": 1, "SIH": 1},
        ],
        psrk_result={"CH3": 4, "CHO": 1, "COO": 1, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)OC(=O)O[Si](C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 5, "CH": 1, "COO": 1, "SIO": 1},
            {"CH3": 5, "CHO": 1, "COO": 1, "SI": 1},
        ],
        psrk_result={"CH3": 5, "CHO": 1, "COO": 1, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH2]OC(=O)OC(C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 4, "C": 1, "COO": 1, "SIH2O": 1},
        psrk_result={"CH3": 4, "C": 1, "COO": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH](C)OC(=O)OC(C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 5, "C": 1, "COO": 1, "SIHO": 1},
        psrk_result={"CH3": 5, "C": 1, "COO": 1, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)(C)OC(=O)O[Si](C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 6, "C": 1, "COO": 1, "SIO": 1},
        psrk_result={"CH3": 6, "C": 1, "COO": 1, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)(C)OC(=O)O[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={},
        psrk_result={},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH]([SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"SIH3": 4, "SIH": 1, "SIH2O": 1, "SIHO": 1},
            {"SIH3": 4, "SIH2": 1, "SIHO": 2},
        ],
        psrk_result={"SIH3": 4, "SIH": 1, "SIH2O": 1, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH]([SiH3])O[SiH2]O[SiH](O[SiH2]O[SiH]([SiH3])[SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"SIH3": 6, "SIH": 1, "SIH2O": 3, "SIHO": 3},
            {"SIH3": 6, "SIH2": 1, "SIH2O": 2, "SIHO": 4},
        ],
        psrk_result={"SIH3": 6, "SIH": 1, "SIH2O": 3, "SIHO": 3},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH]([SiH3])O[SiH2]O[SiH2][SiH](O[SiH2]O[SiH]([SiH3])[SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"SIH3": 6, "SIH": 2, "SIH2O": 4, "SIHO": 2},
            {"SIH3": 6, "SIH2": 1, "SIH": 1, "SIH2O": 3, "SIHO": 3},
            {"SIH3": 6, "SIH2": 2, "SIH2O": 2, "SIHO": 4},
        ],
        psrk_result={"SIH3": 6, "SIH": 2, "SIH2O": 4, "SIHO": 2},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH]([SiH3])O[SiH2]O[SiH]([SiH3])O[SiH2]O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"SIH3": 5, "SIH": 1, "SIH2O": 2, "SIHO": 2},
            {"SIH3": 5, "SIH2": 1, "SIH2O": 1, "SIHO": 3},
        ],
        psrk_result={"SIH3": 5, "SIH": 1, "SIH2O": 2, "SIHO": 2},
        joback_result={},
    ),
    Case(
        identifier="[SiH3]O[SiH2]O[SiH]([SiH3])O[SiH2]O[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={},
        psrk_result={},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH](O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])([SiH3])[SiH3]",  # noqa
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"SIH3": 7, "SI": 1, "SIHO": 1, "SIO": 1},
            {"SIH3": 7, "SIH": 1, "SIO": 2},
        ],
        psrk_result={"SIH3": 7, "SI": 1, "SIHO": 1, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH](O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])(O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3])O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3]",  # noqa
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"SIH3": 13, "SI": 1, "SIHO": 3, "SIO": 3},
            {"SIH3": 13, "SIH": 1, "SIHO": 2, "SIO": 4},
        ],
        psrk_result={"SIH3": 13, "SI": 1, "SIHO": 3, "SIO": 3},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH](O[SiH2][Si]([SiH3])(O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3])O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])([SiH3])[SiH3]",  # noqa
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"SIH3": 13, "SI": 2, "SIH2O": 1, "SIHO": 3, "SIO": 2},
            {"SIH3": 13, "SIH2": 1, "SI": 1, "SIHO": 3, "SIO": 3},
            {"SIH3": 13, "SIH": 1, "SI": 1, "SIH2O": 1, "SIHO": 2, "SIO": 3},
            {"SIH3": 13, "SIH2": 1, "SIH": 1, "SIHO": 2, "SIO": 4},
            {"SIH3": 13, "SIH": 2, "SIH2O": 1, "SIHO": 1, "SIO": 4},
        ],
        psrk_result={"SIH3": 13, "SI": 2, "SIH2O": 1, "SIHO": 3, "SIO": 2},
        joback_result={},
    ),
    Case(
        identifier="[SiH3][SiH](O[Si]([SiH3])([SiH3])[SiH3])O[Si]([SiH3])([SiH3])O[SiH]([SiH3])O[Si]([SiH3])([SiH3])[SiH3]",  # noqa
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"SIH3": 10, "SI": 1, "SIHO": 2, "SIO": 2},
            {"SIH3": 10, "SIH": 1, "SIHO": 1, "SIO": 3},
        ],
        psrk_result={"SIH3": 10, "SI": 1, "SIHO": 2, "SIO": 2},
        joback_result={},
    ),
    Case(
        identifier="[SiH3]OCO[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={},
        psrk_result={},
        joback_result={},
    ),
    Case(
        identifier="COCO[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3O": 1, "CH2O": 1, "SIH3": 1, "SIH2": 1},
            {"CH3": 1, "CH2O": 1, "SIH3": 1, "SIH2O": 1},
            {"CH2": 1, "CH3O": 1, "SIH3": 1, "SIH2O": 1},
        ],
        psrk_result={"CH3O": 1, "CH2O": 1, "SIH3": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3]OCO[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH2O": 1, "SIH3": 2, "SIH2O": 1},
        psrk_result={"CH2O": 1, "SIH3": 2, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH2]OCO[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 1, "CH2O": 1, "SIH3": 1, "SIH2": 1, "SIH2O": 1},
            {"CH3": 1, "CH2": 1, "SIH3": 1, "SIH2O": 2},
        ],
        psrk_result={"CH3": 1, "CH2O": 1, "SIH3": 1, "SIH2": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH](C)OCO[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "CH2O": 1, "SIH3": 1, "SIH2": 1, "SIHO": 1},
            {"CH3": 2, "CH2": 1, "SIH3": 1, "SIH2O": 1, "SIHO": 1},
            {"CH3": 2, "CH2O": 1, "SIH3": 1, "SIH": 1, "SIH2O": 1},
        ],
        psrk_result={"CH3": 2, "CH2O": 1, "SIH3": 1, "SIH": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="C[Si](C)(C)OCO[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CH2O": 1, "SIH3": 1, "SIH2": 1, "SIO": 1},
            {"CH3": 3, "CH2": 1, "SIH3": 1, "SIH2O": 1, "SIO": 1},
            {"CH3": 3, "CH2O": 1, "SIH3": 1, "SI": 1, "SIH2O": 1},
        ],
        psrk_result={"CH3": 3, "CH2O": 1, "SIH3": 1, "SI": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)(C)OCO[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 3, "C": 1, "CH2O": 1, "SIH3": 1, "SIH2O": 1},
        psrk_result={"CH3": 3, "C": 1, "CH2O": 1, "SIH3": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)OCO[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "CH2O": 1, "CHO": 1, "SIH3": 1, "SIH2": 1},
            {"CH3": 2, "CH2": 1, "CHO": 1, "SIH3": 1, "SIH2O": 1},
            {"CH3": 2, "CH": 1, "CH2O": 1, "SIH3": 1, "SIH2O": 1},
        ],
        psrk_result={"CH3": 2, "CH2O": 1, "SIH3": 1, "SIH2": 1, "CHO": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3]OCO[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH2O": 1, "SIH3": 3, "SIHO": 1},
        psrk_result={"CH2O": 1, "SIH3": 3, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="C[SiH](C)OCO[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "CH2O": 1, "SIH3": 2, "SIH": 1, "SIHO": 1},
            {"CH3": 2, "CH2": 1, "SIH3": 2, "SIHO": 2},
        ],
        psrk_result={"CH3": 2, "CH2O": 1, "SIH3": 2, "SIH": 1, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="C[Si](C)(C)OCO[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CH2O": 1, "SIH3": 2, "SIH": 1, "SIO": 1},
            {"CH3": 3, "CH2": 1, "SIH3": 2, "SIHO": 1, "SIO": 1},
            {"CH3": 3, "CH2O": 1, "SIH3": 2, "SI": 1, "SIHO": 1},
        ],
        psrk_result={"CH3": 3, "CH2O": 1, "SIH3": 2, "SI": 1, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)(C)OCO[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 3, "C": 1, "CH2O": 1, "SIH3": 2, "SIHO": 1},
        psrk_result={"CH3": 3, "C": 1, "CH2O": 1, "SIH3": 2, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="COCO[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3O": 1, "CH2O": 1, "SIH3": 2, "SIH": 1},
            {"CH3": 1, "CH2O": 1, "SIH3": 2, "SIHO": 1},
            {"CH2": 1, "CH3O": 1, "SIH3": 2, "SIHO": 1},
        ],
        psrk_result={"CH3O": 1, "CH2O": 1, "SIH3": 2, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)OCO[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "CH2O": 1, "CHO": 1, "SIH3": 2, "SIH": 1},
            {"CH3": 2, "CH2": 1, "CHO": 1, "SIH3": 2, "SIHO": 1},
            {"CH3": 2, "CH": 1, "CH2O": 1, "SIH3": 2, "SIHO": 1},
        ],
        psrk_result={"CH3": 2, "CH2O": 1, "SIH3": 2, "SIH": 1, "CHO": 1},
        joback_result={},
    ),
    Case(
        identifier="CCOCO[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 1, "CH2O": 2, "SIH3": 2, "SIH": 1},
            {"CH3": 1, "CH2": 1, "CH2O": 1, "SIH3": 2, "SIHO": 1},
        ],
        psrk_result={"CH3": 1, "CH2O": 2, "SIH3": 2, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="[SiH3]OCO[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH2O": 1, "SIH3": 4, "SIO": 1},
        psrk_result={"CH2O": 1, "SIH3": 4, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="C[Si](C)(C)OCO[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CH2O": 1, "SIH3": 3, "SI": 1, "SIO": 1},
            {"CH3": 3, "CH2": 1, "SIH3": 3, "SIO": 2},
        ],
        psrk_result={"CH3": 3, "CH2O": 1, "SIH3": 3, "SI": 1, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)(C)OCO[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 3, "C": 1, "CH2O": 1, "SIH3": 3, "SIO": 1},
        psrk_result={"CH3": 3, "C": 1, "CH2O": 1, "SIH3": 3, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="COCO[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3O": 1, "CH2O": 1, "SIH3": 3, "SI": 1},
            {"CH3": 1, "CH2O": 1, "SIH3": 3, "SIO": 1},
            {"CH2": 1, "CH3O": 1, "SIH3": 3, "SIO": 1},
        ],
        psrk_result={"CH3O": 1, "CH2O": 1, "SIH3": 3, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)OCO[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "CH2O": 1, "CHO": 1, "SIH3": 3, "SI": 1},
            {"CH3": 2, "CH2": 1, "CHO": 1, "SIH3": 3, "SIO": 1},
            {"CH3": 2, "CH": 1, "CH2O": 1, "SIH3": 3, "SIO": 1},
        ],
        psrk_result={"CH3": 2, "CH2O": 1, "SIH3": 3, "SI": 1, "CHO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(O[SiH3])O[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 1, "CHO": 1, "SIH3": 2, "SIH2O": 1},
        psrk_result={"CH3": 1, "CHO": 1, "SIH3": 2, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(O[SiH2][SiH3])O[SiH](C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CHO": 1, "SIH3": 1, "SIH": 1, "SIH2O": 1},
            {"CH3": 3, "CH": 1, "SIH3": 1, "SIH2O": 1, "SIHO": 1},
            {"CH3": 3, "CHO": 1, "SIH3": 1, "SIH2": 1, "SIHO": 1},
        ],
        psrk_result={"CH3": 3, "CHO": 1, "SIH3": 1, "SIH": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(O[SiH2][SiH3])OC(C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 4, "C": 1, "CHO": 1, "SIH3": 1, "SIH2O": 1},
        psrk_result={"CH3": 4, "C": 1, "CHO": 1, "SIH3": 1, "SIH2O": 1},
        joback_result={},
    ),
    Case(
        identifier="COC(C)O[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 1, "CH3O": 1, "CHO": 1, "SIH3": 1, "SIH2": 1},
            {"CH3": 2, "CHO": 1, "SIH3": 1, "SIH2O": 1},
            {"CH3": 1, "CH": 1, "CH3O": 1, "SIH3": 1, "SIH2O": 1},
        ],
        psrk_result={"CH3": 1, "CH3O": 1, "CHO": 1, "SIH3": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)OC(C)O[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CHO": 2, "SIH3": 1, "SIH2": 1},
            {"CH3": 3, "CH": 1, "CHO": 1, "SIH3": 1, "SIH2O": 1},
        ],
        psrk_result={"CH3": 3, "CHO": 2, "SIH3": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="CCOC(C)O[SiH2][SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 2, "CH2O": 1, "CHO": 1, "SIH3": 1, "SIH2": 1},
            {"CH3": 2, "CH2": 1, "CHO": 1, "SIH3": 1, "SIH2O": 1},
            {"CH3": 2, "CH": 1, "CH2O": 1, "SIH3": 1, "SIH2O": 1},
        ],
        psrk_result={"CH3": 2, "CH2O": 1, "CHO": 1, "SIH3": 1, "SIH2": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(O[SiH](C)C)O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CHO": 1, "SIH3": 2, "SIH": 1, "SIHO": 1},
            {"CH3": 3, "CH": 1, "SIH3": 2, "SIHO": 2},
        ],
        psrk_result={"CH3": 3, "CHO": 1, "SIH3": 2, "SIH": 1, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(O[SiH3])O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 1, "CHO": 1, "SIH3": 3, "SIHO": 1},
        psrk_result={"CH3": 1, "CHO": 1, "SIH3": 3, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(O[SiH]([SiH3])[SiH3])OC(C)(C)C",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 4, "C": 1, "CHO": 1, "SIH3": 2, "SIHO": 1},
        psrk_result={"CH3": 4, "C": 1, "CHO": 1, "SIH3": 2, "SIHO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)OC(C)O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CHO": 2, "SIH3": 2, "SIH": 1},
            {"CH3": 3, "CH": 1, "CHO": 1, "SIH3": 2, "SIHO": 1},
        ],
        psrk_result={"CH3": 3, "CHO": 2, "SIH3": 2, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="COC(C)O[SiH]([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 1, "CH3O": 1, "CHO": 1, "SIH3": 2, "SIH": 1},
            {"CH3": 2, "CHO": 1, "SIH3": 2, "SIHO": 1},
            {"CH3": 1, "CH": 1, "CH3O": 1, "SIH3": 2, "SIHO": 1},
        ],
        psrk_result={"CH3": 1, "CH3O": 1, "CHO": 1, "SIH3": 2, "SIH": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(O[SiH3])O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 1, "CHO": 1, "SIH3": 4, "SIO": 1},
        psrk_result={"CH3": 1, "CHO": 1, "SIH3": 4, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(OC(C)(C)C)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result={"CH3": 4, "C": 1, "CHO": 1, "SIH3": 3, "SIO": 1},
        psrk_result={"CH3": 4, "C": 1, "CHO": 1, "SIH3": 3, "SIO": 1},
        joback_result={},
    ),
    Case(
        identifier="CC(C)OC(C)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 3, "CHO": 2, "SIH3": 3, "SI": 1},
            {"CH3": 3, "CH": 1, "CHO": 1, "SIH3": 3, "SIO": 1},
        ],
        psrk_result={"CH3": 3, "CHO": 2, "SIH3": 3, "SI": 1},
        joback_result={},
    ),
    Case(
        identifier="COC(C)O[Si]([SiH3])([SiH3])[SiH3]",
        identifier_type="smiles",
        cases_module="silicon",
        r=None,
        q=None,
        unifac_result=[
            {"CH3": 1, "CH3O": 1, "CHO": 1, "SIH3": 3, "SI": 1},
            {"CH3": 2, "CHO": 1, "SIH3": 3, "SIO": 1},
            {"CH3": 1, "CH": 1, "CH3O": 1, "SIH3": 3, "SIO": 1},
        ],
        psrk_result={"CH3": 1, "CH3O": 1, "CHO": 1, "SIH3": 3, "SI": 1},
        joback_result={},
    ),
]
