import os
import sys
from pathlib import Path

import pytest

from ugropy import Groups
from ugropy.writers import to_clapeyron
from ugropy.writers.clapeyron_writers import write_molar_mass


here = Path(__file__).parent
path_db = here / "test_expected_result"

# TODO: Compare differently on MACOS


def test_to_clapeyron():
    with open(path_db / "molarmass.csv", mode="r") as f:
        df_molarmass = f.read()

    with open(path_db / "ogUNIFAC_groups.csv", mode="r") as f:
        df_unifac = f.read()

    with open(path_db / "PSRK_groups.csv", mode="r") as f:
        df_psrk = f.read()

    with open(path_db / "critical.csv", mode="r") as f:
        df_critical = f.read()

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
        df_critical_ugropy = f.read()

    with open(here / "database" / "molarmass.csv", mode="r") as f:
        df_molarmass_ugropy = f.read()

    with open(
        here / "database" / "ogUNIFAC" / "ogUNIFAC_groups.csv", mode="r"
    ) as f:
        df_unifac_ugropy = f.read()

    with open(here / "database" / "PSRK/PSRK_groups.csv", mode="r") as f:
        df_psrk_ugropy = f.read()

    if sys.platform != "darwin":
        assert df_molarmass == df_molarmass_ugropy

    assert df_psrk == df_psrk_ugropy
    assert df_unifac == df_unifac_ugropy

    if sys.platform != "darwin":
        assert df_critical == df_critical_ugropy

    os.remove(here / "database" / "critical.csv")
    os.remove(here / "database" / "molarmass.csv")
    os.remove(here / "database" / "ogUNIFAC" / "ogUNIFAC_groups.csv")
    os.remove(here / "database" / "PSRK" / "PSRK_groups.csv")
    os.rmdir(here / "database" / "ogUNIFAC")
    os.rmdir(here / "database/PSRK")
    os.rmdir(here / "database")


def test_to_clapeyron_batch_name():
    with open(path_db / "molarmass.csv", mode="r") as f:
        df_molarmass = f.read()

    with open(path_db / "ogUNIFAC_groups.csv", mode="r") as f:
        df_unifac = f.read()

    with open(path_db / "PSRK_groups.csv", mode="r") as f:
        df_psrk = f.read()

    with open(path_db / "critical.csv", mode="r") as f:
        df_critical = f.read()

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
        df_critical_ugropy = f.read()

    with open(here / "database" / "otacon_molarmass.csv", mode="r") as f:
        df_molarmass_ugropy = f.read()

    with open(
        here / "database" / "ogUNIFAC" / "otacon_ogUNIFAC_groups.csv", mode="r"
    ) as f:
        df_unifac_ugropy = f.read()

    with open(
        here / "database" / "PSRK" / "otacon_PSRK_groups.csv", mode="r"
    ) as f:
        df_psrk_ugropy = f.read()

    if sys.platform != "darwin":
        assert df_molarmass == df_molarmass_ugropy

    assert df_psrk == df_psrk_ugropy
    assert df_unifac == df_unifac_ugropy

    if sys.platform != "darwin":
        assert df_critical == df_critical_ugropy

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
        df_molarmass_unifac = f.read()

    os.remove(here / "database" / "molarmass.csv")

    # PSRK molar mass
    to_clapeyron(
        ["limonene", "ethanol"],
        psrk_groups=[limonene.psrk.subgroups, ethanol.psrk.subgroups],
        path=(here / "database").resolve(),
    )

    with open(here / "database" / "molarmass.csv", mode="r") as f:
        df_molarmass_psrk = f.read()

    if sys.platform != "darwin":
        assert df_molarmass_unifac == df_molarmass_psrk

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
