"""PSRK groups writer module."""
from io import StringIO
from typing import List

import pandas as pd


def write_psrk(
    path: str,
    batch_name: str,
    molecules_names: List[str],
    psrk_groups: List[dict],
) -> None:
    """Create the DataFrame with the PSRK groups for Clapeyron.jl.

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
    psrk_groups : List[dict], optional
        List of Predictive Soave-Redlich-Kwong groups.

    Returns
    -------
    pd.DataFrame
        DataFrame with the LV-UNIFAC groups for Clapeyron.jl
    """
    data_str = (
        "Clapeyron Database File|\n"
        "PSRK Groups [csvtype = groups,grouptype = PSRK]\n"  # noqa
        "species|groups"
    )

    df = pd.read_csv(StringIO(data_str), sep="|")

    for name, groups in zip(molecules_names, psrk_groups):
        groups_str = "["

        for group in groups.keys():
            groups_str += f'"{group}" => {groups[group]}, '

        groups_str = groups_str[: len(groups_str) - 2]
        groups_str += "]"

        new_row = {
            "Clapeyron Database File": name,
            "Unnamed: 1": f"{groups_str}",
        }

        df.loc[len(df)] = new_row

    df.columns = ["" if "Unnamed" in col else col for col in df.columns]

    # Write .csv
    if batch_name == "":
        write_path = f"{path}/PSRK_groups.csv"
    else:
        write_path = f"{path}/{batch_name}_PSRK_groups.csv"

    with open(write_path, "w", newline="", encoding="utf-8") as file:
        df.to_csv(file, index=False, sep=",")

    with open(write_path, "r") as file:
        unifac = file.read()

    unifac_corrected = unifac.replace('"PSRK', "PSRK")
    unifac_corrected = unifac_corrected.replace('PSRK]"', "PSRK]")

    with open(write_path, "w") as file:
        file.write(unifac_corrected)
