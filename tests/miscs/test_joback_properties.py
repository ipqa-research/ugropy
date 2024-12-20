import numpy as np

import pytest

from pint import Quantity as q

from rdkit import Chem

from ugropy import joback, ureg


@pytest.mark.joback
def test_p_dichlorobenzene():
    mol = joback.get_groups("C1=CC(=CC=C1Cl)Cl", "smiles")

    assert mol.subgroups == {"-Cl": 2, "ring=CH-": 4, "ring=C<": 2}

    assert np.allclose(mol.normal_boiling_point, q(443.4, "K"), atol=1e-2)
    assert mol.normal_boiling_point.units == ureg.kelvin

    assert np.allclose(mol.fusion_temperature, q(256, "K"), atol=1)
    assert mol.fusion_temperature.units == ureg.kelvin

    assert np.allclose(mol.critical_temperature, q(675, "K"), atol=1)
    assert mol.critical_temperature.units == ureg.kelvin

    assert np.allclose(mol.critical_pressure, q(41.5, "bar"), atol=1e-1)
    assert mol.critical_pressure.units == ureg.bar

    assert np.allclose(mol.critical_volume, q(362, "cm^3/mol"), atol=1)
    assert mol.critical_volume.units == ureg.centimeter**3 / ureg.mole

    assert np.allclose(
        mol.ig_enthalpy_formation, q(26.41, "kJ/mol"), atol=1e-2
    )
    assert mol.ig_enthalpy_formation.units == ureg.kilojoule / ureg.mole

    assert np.allclose(mol.ig_gibbs_formation, q(78.56, "kJ/mol"), atol=1e-2)
    assert mol.ig_gibbs_formation.units == ureg.kilojoule / ureg.mole

    assert np.allclose(
        mol.heat_capacity_ideal_gas(298), q(112.3, "J/mol/K"), atol=1
    )
    assert (
        mol.heat_capacity_ideal_gas(298).units
        == ureg.joule / ureg.mole / ureg.kelvin
    )

    assert np.allclose(
        mol.heat_capacity_ideal_gas(400), q(139.2, "J/mol/K"), atol=1
    )
    assert (
        mol.heat_capacity_ideal_gas(400).units
        == ureg.joule / ureg.mole / ureg.kelvin
    )

    assert np.allclose(
        mol.heat_capacity_ideal_gas(800), q(206.8, "J/mol/K"), atol=1
    )
    assert (
        mol.heat_capacity_ideal_gas(800).units
        == ureg.joule / ureg.mole / ureg.kelvin
    )

    assert np.allclose(
        mol.heat_capacity_ideal_gas(1000), q(224.6, "J/mol/K"), atol=1
    )
    assert (
        mol.heat_capacity_ideal_gas(1000).units
        == ureg.joule / ureg.mole / ureg.kelvin
    )

    assert np.allclose(
        mol.vaporization_enthalpy, q(40.66, "kJ/mol"), atol=1e-2
    )
    assert mol.vaporization_enthalpy.units == ureg.kilojoule / ureg.mole

    assert np.allclose(mol.fusion_enthalpy, q(13.3, "kJ/mol"), atol=1e-1)
    assert mol.fusion_enthalpy.units == ureg.kilojoule / ureg.mole

    assert np.allclose(
        mol.viscosity_liquid(333.8), q(7.26e-4, "Pa s"), atol=1e-6
    )
    assert mol.viscosity_liquid(333.8).units == ureg.pascal * ureg.second

    assert np.allclose(
        mol.viscosity_liquid(374.4), q(4.92e-4, "Pa s"), atol=1e-6
    )
    assert mol.viscosity_liquid(374.4).units == ureg.pascal * ureg.second

    assert np.allclose(
        mol.viscosity_liquid(403.1), q(3.91e-4, "Pa s"), atol=1e-6
    )
    assert mol.viscosity_liquid(403.1).units == ureg.pascal * ureg.second

    assert np.allclose(
        mol.viscosity_liquid(423.3), q(3.40e-4, "Pa s"), atol=1e-6
    )
    assert mol.viscosity_liquid(423.3).units == ureg.pascal * ureg.second


@pytest.mark.joback
def test_p_dichlorobenzene_real_nbt():
    mol = joback.get_groups(
        "C1=CC(=CC=C1Cl)Cl", "smiles", normal_boiling_point=447
    )
    assert np.allclose(mol.critical_temperature.magnitude, 681, atol=1)


@pytest.mark.joback
def test_acentric_factor():
    # Perfluorohexane
    mol = joback.get_groups(
        "C(C(C(C(F)(F)F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F",
        "smiles",
        normal_boiling_point=329.8,
    )
    assert np.allclose(mol.acentric_factor.magnitude, 0.525, atol=1e-2)
    assert mol.acentric_factor.units == ureg.dimensionless


@pytest.mark.joback
def test_vapor_pressure():
    mol = joback.get_groups(
        "CC(C)O", identifier_type="smiles", normal_boiling_point=355.4
    )
    assert np.allclose(mol.vapor_pressure(450).magnitude, 15.09, atol=2)
    assert mol.vapor_pressure(450).units == ureg.bar


@pytest.mark.joback
def test_liquid_heat_capacity():
    mol = joback.get_groups("CC(=O)C", identifier_type="smiles")
    assert np.allclose(
        28.0, mol.heat_capacity_liquid(180).magnitude * 0.239, atol=4e-1
    )
    assert (
        mol.heat_capacity_liquid(180).units
        == ureg.joule / ureg.mole / ureg.kelvin
    )

    assert np.allclose(
        28.2, mol.heat_capacity_liquid(209).magnitude * 0.239, atol=1
    )
    assert np.allclose(
        29.8, mol.heat_capacity_liquid(297).magnitude * 0.239, atol=3
    )


@pytest.mark.joback
def test_making_it_explode():
    with pytest.raises(ValueError):
        joback.get_groups("acetone", identifier_type="Messi")


@pytest.mark.joback
def test_rdkit_mol():
    mol1 = joback.get_groups("CCO", identifier_type="smiles")

    chm = Chem.MolFromInchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    mol2 = joback.get_groups(chm, identifier_type="mol")

    assert mol1.subgroups == mol2.subgroups
    assert (
        mol1.experimental_boiling_temperature
        == mol2.experimental_boiling_temperature
    )
    assert mol1.critical_temperature == mol2.critical_temperature
    assert mol1.critical_pressure == mol2.critical_pressure
    assert mol1.critical_volume == mol2.critical_volume
    assert mol1.normal_boiling_point == mol2.normal_boiling_point
    assert mol1.fusion_temperature == mol2.fusion_temperature
    assert mol1.ig_enthalpy_formation == mol2.ig_enthalpy_formation
    assert mol1.ig_gibbs_formation == mol2.ig_gibbs_formation
    assert np.allclose(
        mol1.heat_capacity_ideal_gas_params,
        mol2.heat_capacity_ideal_gas_params,
    )
    assert mol1.fusion_enthalpy == mol2.fusion_enthalpy
    assert mol1.vaporization_enthalpy == mol2.vaporization_enthalpy
    assert mol1.sum_na == mol2.sum_na
    assert mol1.sum_nb == mol2.sum_nb
    assert mol1.molecular_weight == mol2.molecular_weight
    assert mol1.acentric_factor == mol2.acentric_factor
    assert mol1.vapor_pressure_params == mol2.vapor_pressure_params

    assert mol1.heat_capacity_ideal_gas(298) == mol2.heat_capacity_ideal_gas(
        298
    )
    assert mol1.heat_capacity_liquid(298) == mol2.heat_capacity_liquid(298)
    assert mol1.viscosity_liquid(298) == mol2.viscosity_liquid(298)
    assert mol1.vapor_pressure(298) == mol2.vapor_pressure(298)
