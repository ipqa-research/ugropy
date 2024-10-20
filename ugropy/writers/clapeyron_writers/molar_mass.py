"""Molar mass writer module."""

import pathlib
from io import StringIO
from typing import List

import numpy as np

import pandas as pd

from ugropy.core.frag_classes.joback.joback_result import (
    JobackFragmentationResult,
)
from ugropy.models.psrkmod import psrk
from ugropy.models.unifacmod import unifac


def write_molar_mass(
    path: pathlib.Path,
    batch_name: str,
    molecules_names: List[str],
    unifac_groups: List[dict] = [],
    psrk_groups: List[dict] = [],
    joback_objects: List[JobackFragmentationResult] = [],
) -> None:
    """Create the DataFrame with the molecular weights for Clapeyron.jl.

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
    unifac_groups : List[dict], optional
        List of classic liquid-vapor UNIFAC groups, by default [].
    psrk_groups : List[dict], optional
        List of Predictive Soave-Redlich-Kwong groups, by default [].
    joback_objects : List[JobackFragmentationResult], optional
        List of JobackFragmentationResult objects, by default [].

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
    path_molar_mass = pathlib.Path(path)
    # =========================================================================
    # Get molecular weights
    # =========================================================================
    if joback_objects:
        molecular_weigths = [j.molecular_weight for j in joback_objects]
    elif unifac_groups:
        df = unifac.subgroups.copy()
        molecular_weigths = []
        for groups in unifac_groups:
            contribution = df.loc[groups.keys(), "molecular_weight"].to_numpy()
            molecular_weigths.append(
                np.dot(contribution, list(groups.values()))
            )
    elif psrk_groups:
        df = psrk.subgroups.copy()
        molecular_weigths = []
        for groups in psrk_groups:
            contribution = df.loc[groups.keys(), "molecular_weight"].to_numpy()
            molecular_weigths.append(
                np.dot(contribution, list(groups.values()))
            )
    else:
        raise ValueError("Joback, UNIFAC or PSRK groups must be provided.")

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

    df.columns = ["" if "Unnamed" in col else col for col in df.columns]

    if batch_name == "":
        with open(
            path_molar_mass / "molarmass.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as file:
            df.to_csv(file, index=False)

    else:
        with open(
            path_molar_mass / f"{batch_name}_molarmass.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as file:
            df.to_csv(file, index=False)
