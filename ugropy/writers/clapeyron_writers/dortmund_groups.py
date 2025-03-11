"""UNIFAC groups writer module."""

import pathlib
from typing import List


def write_dortmund(
    path: pathlib.Path,
    batch_name: str,
    molecules_names: List[str],
    dortmund_groups: List[dict],
) -> None:
    """Create the DataFrame with the classic LV-UNIFAC groups for Clapeyron.jl.

    Parameters
    ----------
    path : pathlib.Path
        Path to the directory to store de .csv files, by default "./database".
    batch_name : str, optional
        Name of the writing batch. For example, if you name the batch with
        "batch1", the output of the UNIFAC groups will be:
        "batch1_ogUNIFAC_groups.csv". With the default value will be
        "ogUNIFAC_groups.csv", by default "".
    molecules_names : List[str]
        List of names for each chemical to write in the .csv files.
    dortmund_groups : List[dict], optional
        List of classic Dortmund groups.
    """
    lines = [
        "Clapeyron Database File,\n"
        "modified UNIFAC (Dortmund) Groups [csvtype = groups,grouptype = UNIFACDortmund]\n"  # noqa
        "species,groups\n"
    ]

    path_dortmund = path / "Dortmund"

    for name, groups in zip(molecules_names, dortmund_groups):
        groups_str = '"['

        for group in groups.keys():
            groups_str += f'""{group}"" => {groups[group]}, '

        groups_str = groups_str[: len(groups_str) - 2]
        groups_str += ']"\n'

        new_line = [f"{name},{groups_str}"]

        lines.extend(new_line)

    # Create folder for Dortmund groups
    if not path_dortmund.is_dir():
        path_dortmund.mkdir(parents=True)

    # Write .csv
    if batch_name == "":
        write_path = path_dortmund / "UNIFAC_groups.csv"
    else:
        write_path = path_dortmund / f"{batch_name}_UNIFAC_groups.csv"

    with open(write_path, "w", encoding="utf-8", newline="\n") as file:
        file.writelines(lines)
