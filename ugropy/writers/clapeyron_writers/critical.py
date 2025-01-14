"""Critical properties writer module for Clapeyron.jl."""

import pathlib
from io import StringIO
from typing import List

import pandas as pd


def write_critical(
    path: pathlib.Path,
    batch_name: str,
    molecules_names: List[str],
    property_estimator: List = [],
) -> None:
    """Create the DataFrame with the critical properties for Clapeyron.jl.

    Uses the estimated critical properties from any property estimator of
    ugropy.

    Parameters
    ----------
    path : pathlib.Path, optional
        Path to the directory to store de .csv files, by default "./database".
    batch_name : str, optional
        Name of the writing batch. For example, if you name the batch with
        "batch1", the output of the UNIFAC groups will be:
        "batch1_ogUNIFAC_groups.csv". With the default value will be
        "ogUNIFAC_groups.csv", by default "".
    molecules_names : List[str]
        List of names for each chemical to write in the .csv files.
    property_estimator : List, optional
        List of JobackFragmentationResult or AGaniFragmentationResult, by
        default [].
    """
    data_str = (
        "Clapeyron Database File,,,,,\n"
        "Critical Single Parameters,,,,,\n"
        "species,CAS,Tc,Pc,Vc,acentricfactor\n"
    )
    path_critical = pathlib.Path(path)
    # =========================================================================
    # Build dataframe
    # =========================================================================
    df = pd.read_csv(StringIO(data_str))

    for idx, name in enumerate(molecules_names):
        new_row = {
            "Clapeyron Database File": name,
            "Unnamed: 1": "",
            "Unnamed: 2": property_estimator[idx]
            .critical_temperature.to("K")
            .magnitude,
            "Unnamed: 3": property_estimator[idx]
            .critical_pressure.to("Pa")
            .magnitude,
            "Unnamed: 4": property_estimator[idx]
            .critical_volume.to("m^3/mol")
            .magnitude,
            "Unnamed: 5": property_estimator[idx].acentric_factor.magnitude,
        }
        df.loc[len(df)] = new_row

    df.columns = ["" if "Unnamed" in col else col for col in df.columns]

    if batch_name == "":
        with open(
            path_critical / "critical.csv", "w", newline="", encoding="utf-8"
        ) as file:
            df.to_csv(file, index=False)

    else:
        with open(
            path_critical / f"{batch_name}_critical.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as file:
            df.to_csv(file, index=False)
