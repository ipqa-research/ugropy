import pytest

import ugropy as ug


# =============================================================================
# New one molecule groups of PSRK
# =============================================================================
trials_psrk = [
    ("ozone", {"O3": 1}, "name"),
    ("ethylene", {"H2C=CH2": 1}, "name"),
    ("Acetylene", {"CH=-CH": 1}, "name"),
    ("ammonia", {"NH3": 1}, "name"),
    ("carbon monoxide", {"CO": 1}, "name"),
    ("Hydrogen sulfide", {"H2S": 1}, "name"),
    ("nitrogen", {"N2": 1}, "name"),
    ("argon", {"AR": 1}, "name"),
    ("carbon dioxide", {"CO2": 1}, "name"),
    ("methane", {"CH4": 1}, "name"),
    ("oxygen", {"O2": 1}, "name"),
    ("sulfur dioxide", {"SO2": 1}, "name"),
    ("nitric oxide", {"NO": 1}, "name"),
    ("Nitrous oxide", {"N2O": 1}, "name"),
    ("Sulfur hexafluoride", {"SF6": 1}, "name"),
    ("helium", {"HE": 1}, "name"),
    ("neon", {"NE": 1}, "name"),
    ("krypton", {"KR": 1}, "name"),
    ("xenon", {"XE": 1}, "name"),
    ("Hydrogen fluoride", {"HF": 1}, "name"),
    ("hydrogen chloride", {"HCL": 1}, "name"),
    ("Hydrobromic acid", {"HBR": 1}, "name"),
    ("Hydrogen iodide", {"HI": 1}, "name"),
    ("Carbonyl sulfide", {"COS": 1}, "name"),
    ("Ethylene oxide", {"H2COCH2": 1}, "name"),
    ("Fluorine", {"F2": 1}, "name"),
    ("Chlorine", {"CL2": 1}, "name"),
    ("Bromine", {"BR2": 1}, "name"),
    ("Hydrogen cyanide", {"HCN": 1}, "name"),
    ("Nitrogen dioxide", {"NO2": 1}, "name"),
    ("Carbon tetrafluoride", {"CF4": 1}, "name"),
    ("ozone", {"O3": 1}, "name"),
    ("Nitrosyl chloride", {"CLNO": 1}, "name"),
    ("CC(C)(C)N", {"CNH2": 1, "CH3": 3}, "smiles"),
    # TODO
    # ("hydrogen", {"H2": 1}, "name"),
    # ("deuterium", {""}, "name"),
]


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_psrk)
def test_one_group_molecules(identifier, result, identifier_type):
    groups = ug.Groups(identifier, identifier_type, unifac=False)
    assert groups.psrk_groups == result
