from ugropy import abdulelah_gani, dortmund, psrk, unifac


def test_unifac_nonoptimal():
    result = unifac.get_groups(
        "C1OC1C1=CC=CC=C1",
        "smiles",
        search_multiple_solutions=True,
        search_nonoptimal=True,
    )
    assert len(result) == 3

    nonop = [
        {"ACH": 5, "ACCH": 1, "CH2O": 1},
        {"CH2": 1, "ACH": 5, "AC": 1, "CHO": 1},
        {"CH": 1, "ACH": 5, "AC": 1, "CH2O": 1},
    ]

    assert result[0].subgroups in nonop
    assert result[1].subgroups in nonop
    assert result[2].subgroups in nonop


def test_dortmund_nonoptimal():
    result = dortmund.get_groups(
        "C1OC1C1=CC=CC=C1",
        "smiles",
        search_multiple_solutions=True,
        search_nonoptimal=True,
    )
    assert len(result) == 3

    nonop = [
        {"ACH": 5, "AC": 1, "H2COCH": 1},
        {"ACH": 5, "AC": 1, "CHO": 1, "CY-CH2": 1},
        {"ACH": 5, "AC": 1, "CY-CH": 1, "CY-CH2O": 1},
    ]

    assert result[0].subgroups in nonop
    assert result[1].subgroups in nonop
    assert result[2].subgroups in nonop


def test_psrk_nonoptimal():
    result = psrk.get_groups(
        "C1OC1C1=CC=CC=C1",
        "smiles",
        search_multiple_solutions=True,
        search_nonoptimal=True,
    )
    assert len(result) == 4

    nonop = [
        {"ACH": 5, "AC": 1, "H2COCH": 1},
        {"CH2": 1, "ACH": 5, "AC": 1, "CHO": 1},
        {"ACH": 5, "ACCH": 1, "CH2O": 1},
        {"CH": 1, "ACH": 5, "AC": 1, "CH2O": 1},
    ]

    assert result[0].subgroups in nonop
    assert result[1].subgroups in nonop
    assert result[2].subgroups in nonop
    assert result[3].subgroups in nonop


def test_abdulelah_gani_nonoptimal():
    result = abdulelah_gani.get_groups(
        "C1OC1C1=CC=CC=C1",
        "smiles",
        search_multiple_solutions=True,
        search_nonoptimal=True,
    )
    assert len(result) == 2

    assert result[0].primary.subgroups == {
        "aCH": 5,
        "aC except as above": 1,
        "C2H3O": 1,
    }
    assert result[1].primary.subgroups == {
        "aCH": 5,
        "aC except as above": 1,
        "CH2 (cyclic)": 1,
        "CH (cyclic)": 1,
        "O (cyclic)": 1,
    }
