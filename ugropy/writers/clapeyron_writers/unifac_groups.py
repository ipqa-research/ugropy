"""UNIFAC groups writer module."""
import os
from typing import List


def write_unifac(
    path: str,
    batch_name: str,
    molecules_names: List[str],
    unifac_groups: List[dict],
) -> None:
    """Create the DataFrame with the classic LV-UNIFAC groups for Clapeyron.jl.

    Parameters
    ----------
    path : str, optional
        Path to the directory to store de .csv files, by default "./database".
    batch_name : str, optional
        Name of the writing batch. For example, if you name the batch with
        "batch1", the output of the UNIFAC groups will be:
        "batch1_ogUNIFAC_groups.csv". With the default value will be
        "ogUNIFAC_groups.csv", by default "".
    molecules_names : List[str]
        List of names for each chemical to write in the .csv files.
    unifac_groups : List[dict], optional
        List of classic liquid-vapor UNIFAC groups.
    """
    lines = [
        "Clapeyron Database File,\n"
        "original UNIFAC Groups,[csvtype = groups,grouptype = originalUNIFAC]\n"  # noqa
        "species,groups\n"
    ]

    for name, groups in zip(molecules_names, unifac_groups):
        groups_str = '"['

        for group in groups.keys():
            groups_str += f'""{group}"" => {groups[group]}, '

        groups_str = groups_str[: len(groups_str) - 2]
        groups_str += ']"\n'

        new_line = [f"{name},{groups_str}"]

        lines.extend(new_line)

    # Create folder for ogUNIFAC groups
    if not os.path.exists(f"{path}/ogUNIFAC"):
        os.makedirs(f"{path}/ogUNIFAC")

    # Write .csv
    if batch_name == "":
        write_path = f"{path}/ogUNIFAC/ogUNIFAC_groups.csv"
    else:
        write_path = f"{path}/ogUNIFAC/{batch_name}_ogUNIFAC_groups.csv"

    with open(f"{write_path}", "w", encoding="utf-8", newline="\n") as file:
        file.writelines(lines)
