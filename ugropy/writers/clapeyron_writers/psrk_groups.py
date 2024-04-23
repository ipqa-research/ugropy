"""PSRK groups writer module."""

import pathlib
from typing import List


def write_psrk(
    path: pathlib.Path,
    batch_name: str,
    molecules_names: List[str],
    psrk_groups: List[dict],
) -> None:
    """Create the DataFrame with the PSRK groups for Clapeyron.jl.

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
    psrk_groups : List[dict], optional
        List of Predictive Soave-Redlich-Kwong groups.

    Returns
    -------
    pd.DataFrame
        DataFrame with the LV-UNIFAC groups for Clapeyron.jl
    """
    lines = [
        "Clapeyron Database File,\n"
        "PSRK Groups [csvtype = groups,grouptype = PSRK]\n"
        "species,groups\n"
    ]

    path_psrk = path / "PSRK"

    for name, groups in zip(molecules_names, psrk_groups):
        groups_str = '"['

        for group in groups.keys():
            groups_str += f'""{group}"" => {groups[group]}, '

        groups_str = groups_str[: len(groups_str) - 2]
        groups_str += ']"\n'

        new_line = [f"{name},{groups_str}"]

        lines.extend(new_line)

    # Create folder for PSRK groups
    if not path_psrk.is_dir():
        path_psrk.mkdir(parents=True)

    # Write .csv
    if batch_name == "":
        write_path = path_psrk / "PSRK_groups.csv"
    else:
        write_path = path_psrk / f"{batch_name}_PSRK_groups.csv"

    with open(write_path, "w", encoding="utf-8", newline="\n") as file:
        file.writelines(lines)
