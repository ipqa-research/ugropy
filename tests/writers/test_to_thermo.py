from ugropy import get_groups, unifac, writers


def test_to_thermo():
    smiles = ["CCCCCC", "CCC(=O)C"]

    hexa = get_groups(unifac, smiles[0], "smiles")
    thermo_hexa = writers.to_thermo(hexa.subgroups, unifac)

    assert thermo_hexa == {1: 2, 2: 4}

    buta = get_groups(unifac, smiles[1], "smiles")
    thermo_buta = writers.to_thermo(buta.subgroups, unifac)

    assert thermo_buta == {1: 1, 2: 1, 18: 1}
