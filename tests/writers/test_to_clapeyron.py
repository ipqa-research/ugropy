import os
from pathlib import Path

import pandas as pd

import pytest

from ugropy import Groups
from ugropy.writers import to_clapeyron
from ugropy.writers.clapeyron_writers import write_molar_mass


here = Path(__file__).parent
path_db = here / "test_expected_result"


def test_to_clapeyron():
    with open(path_db / "molarmass.csv", mode="r") as f:
        df_molarmass = pd.read_csv(f, sep="|", index_col=None)

    with open(path_db / "ogUNIFAC_groups.csv", mode="r") as f:
        df_unifac = pd.read_csv(f, sep="|", index_col=None)

    with open(path_db / "PSRK_groups.csv", mode="r") as f:
        df_psrk = pd.read_csv(f, sep="|", index_col=None)

    limonene = Groups("CC1=CCC(CC1)C(=C)C", "smiles")
    ethanol = Groups("CCO", "smiles", normal_boiling_temperature=78 + 273.15)

    to_clapeyron(
        ["limonene", "ethanol"],
        [limonene.unifac.subgroups, ethanol.unifac.subgroups],
        [limonene.psrk.subgroups, ethanol.psrk.subgroups],
        [limonene.joback, ethanol.joback],
        (here / "database").resolve(),
    )

    with open(here / "database" / "critical.csv", mode="r") as f:
        df_critical_ugropy = pd.read_csv(f, sep="|", index_col=None)

    with open(here / "database" / "molarmass.csv", mode="r") as f:
        df_molarmass_ugropy = pd.read_csv(f, sep="|", index_col=None)

    with open(
        here / "database" / "ogUNIFAC" / "ogUNIFAC_groups.csv", mode="r"
    ) as f:
        df_unifac_ugropy = pd.read_csv(f, sep="|", index_col=None)

    with open(here / "database" / "PSRK/PSRK_groups.csv", mode="r") as f:
        df_psrk_ugropy = pd.read_csv(f, sep="|", index_col=None)

    assert df_unifac.equals(df_unifac_ugropy)
    assert df_psrk.equals(df_psrk_ugropy)
    assert df_molarmass.equals(df_molarmass_ugropy)

    assert (
        df_critical_ugropy.columns.to_numpy()[0]
        == "Clapeyron Database File,,,,,"
    )
    assert (
        df_critical_ugropy.iloc[0].to_numpy()[0]
        == "Critical Single Parameters,,,,,"
    )
    assert (
        df_critical_ugropy.iloc[1].to_numpy()[0]
        == "species,CAS,Tc,Pc,Vc,acentricfactor"
    )
    assert (
        df_critical_ugropy.iloc[2].to_numpy()[0]
        == "limonene,,657.4486692170663,2755561.066677689,0.0004965,0.3219127551737521"  # noqa
    )
    assert (
        df_critical_ugropy.iloc[3].to_numpy()[0]
        == "ethanol,,519.5440517502864,5756641.437226128,0.00016649999999999998,0.5610678012451401"  # noqa
    )

    os.remove(here / "database" / "critical.csv")
    os.remove(here / "database" / "molarmass.csv")
    os.remove(here / "database" / "ogUNIFAC" / "ogUNIFAC_groups.csv")
    os.remove(here / "database" / "PSRK" / "PSRK_groups.csv")
    os.rmdir(here / "database" / "ogUNIFAC")
    os.rmdir(here / "database/PSRK")
    os.rmdir(here / "database")


def test_to_clapeyron_batch_name():
    with open(path_db / "molarmass.csv", mode="r") as f:
        df_molarmass = pd.read_csv(f, sep="|", index_col=None)

    with open(path_db / "ogUNIFAC_groups.csv", mode="r") as f:
        df_unifac = pd.read_csv(f, sep="|", index_col=None)

    with open(path_db / "PSRK_groups.csv", mode="r") as f:
        df_psrk = pd.read_csv(f, sep="|", index_col=None)

    limonene = Groups("CC1=CCC(CC1)C(=C)C", "smiles")
    ethanol = Groups("CCO", "smiles", normal_boiling_temperature=78 + 273.15)

    to_clapeyron(
        ["limonene", "ethanol"],
        [limonene.unifac.subgroups, ethanol.unifac.subgroups],
        [limonene.psrk.subgroups, ethanol.psrk.subgroups],
        [limonene.joback, ethanol.joback],
        (here / "database").resolve(),
        "otacon",
    )

    with open(here / "database" / "otacon_critical.csv", mode="r") as f:
        df_critical_ugropy = pd.read_csv(f, sep="|", index_col=None)

    with open(here / "database" / "otacon_molarmass.csv", mode="r") as f:
        df_molarmass_ugropy = pd.read_csv(f, sep="|", index_col=None)

    with open(
        here / "database" / "ogUNIFAC" / "otacon_ogUNIFAC_groups.csv", mode="r"
    ) as f:
        df_unifac_ugropy = pd.read_csv(f, sep="|", index_col=None)

    with open(
        here / "database" / "PSRK" / "otacon_PSRK_groups.csv", mode="r"
    ) as f:
        df_psrk_ugropy = pd.read_csv(f, sep="|", index_col=None)

    assert df_unifac.equals(df_unifac_ugropy)
    assert df_psrk.equals(df_psrk_ugropy)
    assert df_molarmass.equals(df_molarmass_ugropy)

    assert (
        df_critical_ugropy.columns.to_numpy()[0]
        == "Clapeyron Database File,,,,,"
    )
    assert (
        df_critical_ugropy.iloc[0].to_numpy()[0]
        == "Critical Single Parameters,,,,,"
    )
    assert (
        df_critical_ugropy.iloc[1].to_numpy()[0]
        == "species,CAS,Tc,Pc,Vc,acentricfactor"
    )
    assert (
        df_critical_ugropy.iloc[2].to_numpy()[0]
        == "limonene,,657.4486692170663,2755561.066677689,0.0004965,0.3219127551737521"  # noqa
    )
    assert (
        df_critical_ugropy.iloc[3].to_numpy()[0]
        == "ethanol,,519.5440517502864,5756641.437226128,0.00016649999999999998,0.5610678012451401"  # noqa
    )

    os.remove(here / "database" / "otacon_critical.csv")
    os.remove(here / "database" / "otacon_molarmass.csv")
    os.remove(here / "database" / "ogUNIFAC" / "otacon_ogUNIFAC_groups.csv")
    os.remove(here / "database" / "PSRK" / "otacon_PSRK_groups.csv")
    os.rmdir(here / "database" / "ogUNIFAC")
    os.rmdir(here / "database" / "PSRK")
    os.rmdir(here / "database")


def test_molar_mass_csv():

    limonene = Groups("CC1=CCC(CC1)C(=C)C", "smiles")
    ethanol = Groups("CCO", "smiles", normal_boiling_temperature=78 + 273.15)

    # UNIFAC molar mass
    to_clapeyron(
        ["limonene", "ethanol"],
        unifac_groups=[limonene.unifac.subgroups, ethanol.unifac.subgroups],
        psrk_groups=[limonene.psrk.subgroups, ethanol.psrk.subgroups],
        path=(here / "database").resolve(),
    )

    with open(here / "database" / "molarmass.csv", mode="r") as f:
        df_molarmass_unifac = pd.read_csv(f, sep="|", index_col=None)

    os.remove(here / "database" / "molarmass.csv")

    # PSRK molar mass
    to_clapeyron(
        ["limonene", "ethanol"],
        psrk_groups=[limonene.psrk.subgroups, ethanol.psrk.subgroups],
        path=(here / "database").resolve(),
    )

    with open(here / "database" / "molarmass.csv", mode="r") as f:
        df_molarmass_psrk = pd.read_csv(f, sep="|", index_col=None)

    assert (
        df_molarmass_unifac.columns.to_numpy()[0]
        == "Clapeyron Database File,,"
    )
    assert (
        df_molarmass_unifac.iloc[0].to_numpy()[0]
        == "Molar Mases Single Params,,"
    )
    assert df_molarmass_unifac.iloc[1].to_numpy()[0] == "species,CAS,Mw"
    assert df_molarmass_unifac.iloc[2].to_numpy()[0] == "limonene,,136.238"
    assert df_molarmass_unifac.iloc[3].to_numpy()[0] == "ethanol,,46.069"
    assert df_molarmass_unifac.equals(df_molarmass_psrk)

    os.remove(here / "database" / "molarmass.csv")
    os.remove(here / "database" / "ogUNIFAC" / "ogUNIFAC_groups.csv")
    os.remove(here / "database" / "PSRK" / "PSRK_groups.csv")
    os.rmdir(here / "database" / "ogUNIFAC")
    os.rmdir(here / "database" / "PSRK")
    os.rmdir(here / "database")

    # Making it fail
    with pytest.raises(ValueError):
        write_molar_mass(
            path=(here / "database").resolve(),
            batch_name=None,
            molecules_names=["limonene", "ethanol"],
        )


def test_making_it_explode():
    with pytest.raises(ValueError):
        to_clapeyron(molecules_names=[])

    with pytest.raises(ValueError):
        to_clapeyron(
            molecules_names=["limonene", "ethanol"],
            unifac_groups=[{"CH3": 1, "CH2": 1, "OH": 1}],
        )

    with pytest.raises(ValueError):
        to_clapeyron(
            molecules_names=["limonene", "ethanol"],
            psrk_groups=[{"CH3": 1, "CH2": 1, "OH": 1}],
        )

    with pytest.raises(ValueError):
        to_clapeyron(
            molecules_names=["limonene", "ethanol"], joback_objects=["hello"]
        )
