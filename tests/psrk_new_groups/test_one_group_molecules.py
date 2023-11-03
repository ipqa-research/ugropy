import pytest

import ugropy as ug


# =============================================================================
# New one molecule groups of PSRK
# =============================================================================
trials_psrk = [
    # ozone
    ("[O-][O+]=O", {"O3": 1}, "smiles"),
    # ethylene
    ("C=C", {"H2C=CH2": 1}, "smiles"),
    # Acetylene
    ("C#C", {"CH=-CH": 1}, "smiles"),
    # ammonia
    ("N", {"NH3": 1}, "smiles"),
    # carbon monoxide
    ("[C-]#[O+]", {"CO": 1}, "smiles"),
    # Hydrogen sulfide
    ("S", {"H2S": 1}, "smiles"),
    # nitrogen
    ("N#N", {"N2": 1}, "smiles"),
    # argon
    ("[Ar]", {"AR": 1}, "smiles"),
    # carbon dioxide
    ("C(=O)=O", {"CO2": 1}, "smiles"),
    # methane
    ("C", {"CH4": 1}, "smiles"),
    # oxygen
    ("O=O", {"O2": 1}, "smiles"),
    # sulfur dioxide
    ("O=S=O", {"SO2": 1}, "smiles"),
    # nitric oxide
    ("[N]=O", {"NO": 1}, "smiles"),
    # Nitrous oxide
    ("[N-]=[N+]=O", {"N2O": 1}, "smiles"),
    # Sulfur hexafluoride
    ("FS(F)(F)(F)(F)F", {"SF6": 1}, "smiles"),
    # helium
    ("[He]", {"HE": 1}, "smiles"),
    # neon
    ("[Ne]", {"NE": 1}, "smiles"),
    # krypton
    ("[Kr]", {"KR": 1}, "smiles"),
    # xenon
    ("[Xe]", {"XE": 1}, "smiles"),
    # Hydrogen fluoride
    ("F", {"HF": 1}, "smiles"),
    # hydrogen chloride
    ("Cl", {"HCL": 1}, "smiles"),
    # Hydrobromic acid
    ("Br", {"HBR": 1}, "smiles"),
    # Hydrogen iodide
    ("I", {"HI": 1}, "smiles"),
    # Carbonyl sulfide
    ("C(=O)=S", {"COS": 1}, "smiles"),
    # Ethylene oxide
    ("C1CO1", {"H2COCH2": 1}, "smiles"),
    # Fluorine
    ("FF", {"F2": 1}, "smiles"),
    # Chlorine
    ("ClCl", {"CL2": 1}, "smiles"),
    # Bromine
    ("BrBr", {"BR2": 1}, "smiles"),
    # Hydrogen cyanide
    ("C#N", {"HCN": 1}, "smiles"),
    # Nitrogen dioxide
    ("N(=O)[O]", {"NO2": 1}, "smiles"),
    # Carbon tetrafluoride
    ("C(F)(F)(F)F", {"CF4": 1}, "smiles"),
    # Ozone
    ("[O-][O+]=O", {"O3": 1}, "smiles"),
    # Nitrosyl chloride
    ("N(=O)Cl", {"CLNO": 1}, "smiles"),
    # tert-Butylamine
    ("CC(C)(C)N", {"CNH2": 1, "CH3": 3}, "smiles"),
    # TODO
    # ("hydrogen", {"H2": 1}, "name"),
    # ("deuterium", {""}, "name"),
]


@pytest.mark.PSRK
@pytest.mark.parametrize("identifier, result, identifier_type", trials_psrk)
def test_one_group_molecules(identifier, result, identifier_type):
    assert ug.get_psrk_groups(identifier, identifier_type) == result


# =============================================================================
# UNIFAC
# =============================================================================
trials_unifac = [
    # Ethylene oxide
    ("C1CO1", {"CH2O": 1, "CH2": 1}, "smiles"),
]


@pytest.mark.UNIFAC
@pytest.mark.parametrize("identifier, result, identifier_type", trials_unifac)
def test_one_group_molecules_unifac(identifier, result, identifier_type):
    assert ug.get_unifac_groups(identifier, identifier_type) == result
