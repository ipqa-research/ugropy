"""to_clapeyron module."""

import pathlib
from typing import List

from ugropy.writers.clapeyron_writers.critical import write_critical
from ugropy.writers.clapeyron_writers.dortmund_groups import write_dortmund
from ugropy.writers.clapeyron_writers.molar_mass import write_molar_mass
from ugropy.writers.clapeyron_writers.psrk_groups import write_psrk
from ugropy.writers.clapeyron_writers.unifac_groups import write_unifac


def to_clapeyron(
    molecules_names: List[str],
    unifac_groups: List[dict] = [],
    psrk_groups: List[dict] = [],
    dortmund_groups: List[dict] = [],
    property_estimator: List = [],
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
    dortmund_groups : List[dict], optional
        List of Dortmund UNIFAC groups, by default [].
    property_estimator : List, optional
        List of JobackFragmentationResult or AGaniFragmentationResult, by
        default [].
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

    if dortmund_groups and len(dortmund_groups) != len(molecules_names):
        raise ValueError(
            "Dortmund groups list must have the same amount of elements than"
            "the molecules name list."
        )

    if property_estimator and len(property_estimator) != len(molecules_names):
        raise ValueError(
            "Property estimators result objects list must have the same amount"
            "of elements than the molecules name list."
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
        dortmund_groups,
        property_estimator,
    )

    # LV-UNIFAC
    if unifac_groups:
        write_unifac(path_pathlib, batch_name, molecules_names, unifac_groups)

    # PSRK
    if psrk_groups:
        write_psrk(path_pathlib, batch_name, molecules_names, psrk_groups)

    # Dortmund
    if dortmund_groups:
        write_dortmund(
            path_pathlib, batch_name, molecules_names, dortmund_groups
        )

    # Critical
    if property_estimator:
        write_critical(
            path_pathlib, batch_name, molecules_names, property_estimator
        )
