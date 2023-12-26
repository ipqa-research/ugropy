"""Joback critical properties writer module."""
from io import StringIO
from typing import List

import pandas as pd

from ugropy.joback import Joback


def write_critical(
    path: str,
    batch_name: str,
    molecules_names: List[str],
    joback_objects: List[Joback] = [],
) -> None:
    """Create the DataFrame with the critical properties for Clapeyron.jl.

    Uses the Joback to estimate the critical properties of the molecules.

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
    joback_objects : List[Joback], optional
        List of ugropy.Joback objects, by default [].

    Returns
    -------
    pd.DataFrame
        DataFrame with the molecular weights for Clapeyron.jl
    """
    data_str = (
        "Clapeyron Database File,,,,,\n"
        "Critical Single Parameters,,,,,\n"
        "species,CAS,Tc,Pc,Vc,acentricfactor\n"
    )
    # =========================================================================
    # Build dataframe
    # =========================================================================
    df = pd.read_csv(StringIO(data_str))

    for idx, name in enumerate(molecules_names):
        new_row = {
            "Clapeyron Database File": name,
            "Unnamed: 1": "",
            "Unnamed: 2": joback_objects[idx].critical_temperature,
            "Unnamed: 3": joback_objects[idx].critical_pressure * 1e5,
            "Unnamed: 4": joback_objects[idx].critical_volume * 1e-6,
            "Unnamed: 5": joback_objects[idx].acentric_factor,
        }
        df.loc[len(df)] = new_row

    df.columns = ["" if "Unnamed" in col else col for col in df.columns]

    if batch_name == "":
        with open(
            f"{path}/critical.csv", "w", newline="", encoding="utf-8"
        ) as file:
            df.to_csv(file, index=False)

    else:
        with open(
            f"{path}/{batch_name}_critical.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as file:
            df.to_csv(file, index=False)
