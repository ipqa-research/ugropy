import pytest

from ugropy import get_groups, joback


# =============================================================================
# -F, -Cl, -Br, -I
# =============================================================================
# Joback
trials = [
    (
        "ClC(I)=C(Br)C=CC=C",
        {"=CH2": 1, "=CH-": 3, "=C<": 2, "-I": 1, "-Br": 1, "-Cl": 1},
        "smiles",
    ),
    # trichloroethylene
    ("C(=C(Cl)Cl)Cl", {"=CH-": 1, "=C<": 1, "-Cl": 3}, "smiles"),
    # BTI
    ("FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F", {}, "smiles"),
    # 1-chlorobutane
    ("CCCCCl", {"-CH3": 1, "-CH2-": 3, "-Cl": 1}, "smiles"),
    # 2-chloropropane
    ("CC(C)Cl", {"-CH3": 2, ">CH-": 1, "-Cl": 1}, "smiles"),
    # 2-chloro-2-methylpropane
    ("CC(C)(C)Cl", {"-CH3": 3, ">C<": 1, "-Cl": 1}, "smiles"),
    ("OC(Cl)Cl", {">CH-": 1, "-Cl": 2, "-OH (alcohol)": 1}, "smiles"),
    # dichloro methane
    ("C(Cl)Cl", {"-CH2-": 1, "-Cl": 2}, "smiles"),
    # 1,1-dichloroethane
    ("CC(Cl)Cl", {"-CH3": 1, ">CH-": 1, "-Cl": 2}, "smiles"),
    # 2,2-dichloropropane
    ("CC(C)(Cl)Cl", {"-CH3": 2, ">C<": 1, "-Cl": 2}, "smiles"),
    # chloroform
    ("C(Cl)(Cl)Cl", {">CH-": 1, "-Cl": 3}, "smiles"),
    # trichloroethane
    ("CC(Cl)(Cl)Cl", {"-CH3": 1, ">C<": 1, "-Cl": 3}, "smiles"),
    # tetrachloromethane
    ("C(Cl)(Cl)(Cl)Cl", {">C<": 1, "-Cl": 4}, "smiles"),
    # chlorobenzene
    ("C1=CC=C(C=C1)Cl", {"ring=CH-": 5, "-Cl": 1, "ring=C<": 1}, "smiles"),
    # 1-iodoethane
    ("CCI", {"-CH3": 1, "-CH2-": 1, "-I": 1}, "smiles"),
    # Iodobenzene
    ("C1=CC=C(C=C1)I", {"ring=CH-": 5, "ring=C<": 1, "-I": 1}, "smiles"),
    # 1-bromoethane
    ("CCBr", {"-CH3": 1, "-CH2-": 1, "-Br": 1}, "smiles"),
    # Bromobenzene
    ("C1=CC=C(C=C1)Br", {"ring=CH-": 5, "ring=C<": 1, "-Br": 1}, "smiles"),
    # hexafluorobenzene
    ("C1(=C(C(=C(C(=C1F)F)F)F)F)F", {"ring=C<": 6, "-F": 6}, "smiles"),
    (
        "FC1=CC=NC=C1",
        {"ring=CH-": 4, "ring=C<": 1, "-N= (ring)": 1, "-F": 1},
        "smiles",
    ),
    (
        "OC(F)(Br)I",
        {">C<": 1, "-F": 1, "-Br": 1, "-I": 1, "-OH (alcohol)": 1},
        "smiles",
    ),
    ("OC(O)(F)F", {">C<": 1, "-F": 2, "-OH (alcohol)": 2}, "smiles"),
    ("OC(F)(F)F", {">C<": 1, "-F": 3, "-OH (alcohol)": 1}, "smiles"),
    # Perfluorohexane
    (
        "C(C(C(C(F)(F)F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F",
        {">C<": 6, "-F": 14},
        "smiles",
    ),
    # PerfluoromethylcyClohexane
    (
        "C1(C(C(C(C(C1(F)F)(F)F)(F)F)(F)F)(F)F)(C(F)(F)F)F",
        {">C<": 1, "ring>C<": 6, "-F": 14},
        "smiles",
    ),
    ("FC(F)F", {"-F": 3, ">CH-": 1}, "smiles"),
    ("FCF", {"-F": 2, "-CH2-": 1}, "smiles"),
    ("CF", {"-F": 1, "-CH3": 1}, "smiles"),
    # Trichlorofluoromethane
    ("C(F)(Cl)(Cl)Cl", {">C<": 1, "-Cl": 3, "-F": 1}, "smiles"),
    # Tetrachloro-1,2-difluoroethane
    ("C(C(F)(Cl)Cl)(F)(Cl)Cl", {">C<": 2, "-Cl": 4, "-F": 2}, "smiles"),
    # Dichlorofluoromethane
    ("C(F)(Cl)Cl", {">CH-": 1, "-Cl": 2, "-F": 1}, "smiles"),
    # 1-Chloro-1,2,2,2-tetrafluoroethane
    ("C(C(F)(F)F)(F)Cl", {">CH-": 1, ">C<": 1, "-Cl": 1, "-F": 4}, "smiles"),
    # 1,2-Dichlorotetrafluoroethane
    ("C(C(F)(F)Cl)(F)(F)Cl", {">C<": 2, "-Cl": 2, "-F": 4}, "smiles"),
    # Chlorodifluoromethane
    ("C(F)(F)Cl", {">CH-": 1, "-Cl": 1, "-F": 2}, "smiles"),
    # Chlorotrifluoromethane
    ("C(F)(F)(F)Cl", {">C<": 1, "-Cl": 1, "-F": 3}, "smiles"),
    # Dichlorodifluoromethane
    ("C(F)(F)(Cl)Cl", {">C<": 1, "-Cl": 2, "-F": 2}, "smiles"),
]


@pytest.mark.Joback
@pytest.mark.parametrize("identifier, result, identifier_type", trials)
def test_joback_halogens(identifier, result, identifier_type):
    assert get_groups(joback, identifier, identifier_type).subgroups == result
