"""to_clapeyron module."""

import pathlib
from typing import List

from ugropy.properties.joback_properties import JobackProperties

from .clapeyron_writers import (
    write_critical,
    write_molar_mass,
    write_psrk,
    write_unifac,
)


def to_clapeyron(
    molecules_names: List[str],
    unifac_groups: List[dict] = [],
    psrk_groups: List[dict] = [],
    joback_objects: List[JobackProperties] = [],
    path: str = "database",
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
    joback_objects : List[JobackProperties], optional
        List of ugropy.properties.JobackProperties objects, by default [].
    path : str, optional
        Path to the directory to store de .csv files, by default "./database".
    batch_name : str, optional
        Name of the writing batch. For example, if you name the batch with
        "batch1", the output of the UNIFAC groups will be:
        "batch1_ogUNIFAC_groups.csv". With the default value will be
        "ogUNIFAC_groups.csv", by default "".
    """
    # Use pathlib's Path internally
    path_pathlib = pathlib.Path(path)

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

    # Create dir if not created
    if not path_pathlib.is_dir():
        path_pathlib.mkdir(parents=True)

    # Molar mass
    write_molar_mass(
        path_pathlib,
        batch_name,
        molecules_names,
        unifac_groups,
        psrk_groups,
        joback_objects,
    )

    # LV-UNIFAC
    if unifac_groups:
        write_unifac(path_pathlib, batch_name, molecules_names, unifac_groups)

    # PSRK
    if psrk_groups:
        write_psrk(path_pathlib, batch_name, molecules_names, psrk_groups)

    # Critical
    if joback_objects:
        write_critical(
            path_pathlib, batch_name, molecules_names, joback_objects
        )
