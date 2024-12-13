from ugropy import unifac, writers


def test_to_thermo():
    smiles = ["CCCCCC", "CCC(=O)C"]

    hexa = unifac.get_groups(smiles[0], "smiles")
    thermo_hexa = writers.to_thermo(hexa.subgroups, unifac)

    assert thermo_hexa == {1: 2, 2: 4}

    buta = unifac.get_groups(
        smiles[1], "smiles", search_multiple_solutions=True
    )
    thermo_buta = [writers.to_thermo(b.subgroups, unifac) for b in buta]

    sols = [{1: 1, 2: 1, 18: 1}, {1: 2, 19: 1}]

    assert thermo_buta[0] in sols
    assert thermo_buta[1] in sols
