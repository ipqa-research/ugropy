from ugropy import dortmund, psrk, unifac, writers


def test_to_yaeos_dortmund():
    identifiers = ["CC", "CCO", "C1CCCCC1", "O"]

    groups = [dortmund.get_groups(iden, "smiles") for iden in identifiers]

    fortran_code = writers.to_yaeos(groups)

    expected = (
        "use yaeos__models_ge_group_contribution_unifac, only: Groups\n"
        "\n"
        "type(Groups) :: molecules(4)\n"
        "\n"
        "molecules(1)%groups_ids = [1]\n"
        "molecules(1)%number_of_groups = [2]\n"
        "\n"
        "molecules(2)%groups_ids = [1, 2, 14]\n"
        "molecules(2)%number_of_groups = [1, 1, 1]\n"
        "\n"
        "molecules(3)%groups_ids = [78]\n"
        "molecules(3)%number_of_groups = [6]\n"
        "\n"
        "molecules(4)%groups_ids = [16]\n"
        "molecules(4)%number_of_groups = [1]\n"
        "\n"
    )

    assert fortran_code == expected


def test_to_yaeos_psrk():
    identifiers = ["CC", "CCO", "C1CCCCC1", "O=O"]

    groups = [psrk.get_groups(iden, "smiles") for iden in identifiers]

    fortran_code = writers.to_yaeos(groups)

    expected = (
        "use yaeos__models_ge_group_contribution_unifac, only: Groups\n"
        "\n"
        "type(Groups) :: molecules(4)\n"
        "\n"
        "molecules(1)%groups_ids = [1]\n"
        "molecules(1)%number_of_groups = [2]\n"
        "\n"
        "molecules(2)%groups_ids = [1, 2, 14]\n"
        "molecules(2)%number_of_groups = [1, 1, 1]\n"
        "\n"
        "molecules(3)%groups_ids = [2]\n"
        "molecules(3)%number_of_groups = [6]\n"
        "\n"
        "molecules(4)%groups_ids = [119]\n"
        "molecules(4)%number_of_groups = [1]\n"
        "\n"
    )

    assert fortran_code == expected


def test_to_yaeos_unifac():
    identifiers = ["CC", "CCO", "c1ccccc1C", "O"]

    groups = [unifac.get_groups(iden, "smiles") for iden in identifiers]

    fortran_code = writers.to_yaeos(groups)

    expected = (
        "use yaeos__models_ge_group_contribution_unifac, only: Groups\n"
        "\n"
        "type(Groups) :: molecules(4)\n"
        "\n"
        "molecules(1)%groups_ids = [1]\n"
        "molecules(1)%number_of_groups = [2]\n"
        "\n"
        "molecules(2)%groups_ids = [1, 2, 14]\n"
        "molecules(2)%number_of_groups = [1, 1, 1]\n"
        "\n"
        "molecules(3)%groups_ids = [9, 11]\n"
        "molecules(3)%number_of_groups = [5, 1]\n"
        "\n"
        "molecules(4)%groups_ids = [16]\n"
        "molecules(4)%number_of_groups = [1]\n"
        "\n"
    )

    assert fortran_code == expected
