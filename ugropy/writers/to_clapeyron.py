"""to_clapeyron module."""
from io import StringIO
from typing import List

import numpy as np
from numpy.typing import NDArray

from ugropy import Joback

import pandas as pd


def to_clapeyron(
    molecules_names: List[str],
    unifac_groups: List[dict] = [],
    psrk_groups: List[dict] = [],
    joback_objects: List[Joback] = [],
    path: str = "./database",
    batch_name: str = "",
) -> None:
    """Write the .csv input files for Clapeyron.jl.

    The provided lists must have the same length. If one of the model lists is
    left empty, that model will not be writen in the .csv.

    Parameters
    ----------
    molecules_names : List[str]
        List of names for each chemical to write in the .csv files.
    unifac_groups : List[dict], optional
        List of classic liquid-vapor UNIFAC groups, by default [].
    psrk_groups : List[dict], optional
        List of Predictive Soave-Redlich-Kwong groups, by default [].
    joback_objects : List[Joback], optional
        List of ugropy.Joback objects, by default [].
    path : str, optional
        Path to the directory to store de .csv files, by default "./database".
    batch_name : str, optional
        Name of the writing batch. For example, if you name the batch with
        "batch1", the output of the UNIFAC groups will be:
        "batch1_ogUNIFAC_groups.csv". With the default value will be
        "ogUNIFAC_groups.csv", by default "".
    """
    # Check if all list have correct data:
    if len(molecules_names) == 0:
        raise ValueError("No names provided for the molecules.")

    if unifac_groups and len(unifac_groups) != len(molecules_names):
        raise ValueError(
            "UNIFAC groups list must have the same amount of elements than"
            "the molecules name list."
        )

    if psrk_groups and len(psrk_groups) != len(molecules_names):
        raise ValueError(
            "PSRK groups list must have the same amount of elements than"
            "the molecules name list."
        )

    if joback_objects and len(joback_objects) != len(molecules_names):
        raise ValueError(
            "Joback objects list must have the same amount of elements than"
            "the molecules name list."
        )

    # Molar mass
    molarmass_df = _get_molar_mass_df(
        molecules_names, unifac_groups, psrk_groups, joback_objects
    )

    molarmass_df.to_csv(f"{path}/molarmass.csv", index=False)


def _get_molar_mass_df(
    molecules_names: List[str],
    unifac_groups: List[dict] = [],
    psrk_groups: List[dict] = [],
    joback_objects: List[Joback] = [],
) -> pd.DataFrame:
    """Create the DataFrame with the molecular weights for Clapeyron.jl.

    Parameters
    ----------
    molecules_names : List[str]
        List of names for each chemical to write in the .csv files.
    unifac_groups : List[dict], optional
        List of classic liquid-vapor UNIFAC groups, by default [].
    psrk_groups : List[dict], optional
        List of Predictive Soave-Redlich-Kwong groups, by default [].
    joback_objects : List[Joback], optional
        List of ugropy.Joback objects, by default [].

    Returns
    -------
    pd.DataFrame
        DataFrame with the molecular weights for Clapeyron.jl
    """
    data_str = (
        "Clapeyron Database File,,\n"
        "Molar Mases Single Params,,\n"
        "species,CAS,Mw\n"
    )
    # =========================================================================
    # Get molecular weights
    # =========================================================================
    if joback_objects:
        molecular_weigths = [j.molecular_weight for j in joback_objects]

    # =========================================================================
    # Build dataframe
    # =========================================================================
    df = pd.read_csv(StringIO(data_str))

    for idx, name in enumerate(molecules_names):
        new_row = {
            "Clapeyron Database File": name,
            "Unnamed: 1": "",
            "Unnamed: 2": molecular_weigths[idx],
        }
        df.loc[len(df)] = new_row

    # explicar linea
    df.columns = ["" if "Unnamed" in col else col for col in df.columns]

    return df