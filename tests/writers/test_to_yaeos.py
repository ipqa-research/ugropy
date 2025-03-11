from ugropy import dortmund, psrk, unifac, writers


def test_to_yaeos_dortmund():
    names = ["ethane", "ethanol", "cyclohexane", "water"]

    groups = [dortmund.get_groups(name).subgroups for name in names]

    fortran_code = writers.to_yaeos(groups, dortmund)

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
    names = ["ethane", "ethanol", "cyclohexane", "oxygen"]

    groups = [psrk.get_groups(name).subgroups for name in names]

    fortran_code = writers.to_yaeos(groups, psrk)

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
    names = ["ethane", "ethanol", "toluene", "water"]

    groups = [unifac.get_groups(name).subgroups for name in names]

    fortran_code = writers.to_yaeos(groups, unifac)

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
