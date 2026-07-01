from pathlib import Path

import ugropy as ug

expected = Path(__file__).parent / "expected.txt"


def test_expected_svg():
    with open(expected, "r") as f:
        expected_content = f.read()

    result = ug.unifac.get_groups("CCCC", "smiles")

    svg = result.get_solution_svg()

    assert svg == expected_content
