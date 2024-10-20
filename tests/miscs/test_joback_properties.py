import numpy as np

import pytest

from rdkit import Chem

from ugropy import joback


@pytest.mark.joback
def test_p_dichlorobenzene():
    mol = joback.get_groups("C1=CC(=CC=C1Cl)Cl", "smiles")
    assert mol.subgroups == {"-Cl": 2, "ring=CH-": 4, "ring=C<": 2}
    assert np.allclose(mol.normal_boiling_point, 443.4, atol=1e-2)
    assert np.allclose(mol.fusion_temperature, 256, atol=1)
    assert np.allclose(mol.critical_temperature, 675, atol=1)
    assert np.allclose(mol.critical_pressure, 41.5, atol=1e-1)
    assert np.allclose(mol.critical_volume, 362, atol=1)
    assert np.allclose(mol.h_formation, 26.41, atol=1e-2)
    assert np.allclose(mol.g_formation, 78.56, atol=1e-2)
    assert np.allclose(mol.heat_capacity_ideal_gas(298), 112.3, atol=1)
    assert np.allclose(mol.heat_capacity_ideal_gas(400), 139.2, atol=1)
    assert np.allclose(mol.heat_capacity_ideal_gas(800), 206.8, atol=1)
    assert np.allclose(mol.heat_capacity_ideal_gas(1000), 224.6, atol=1)
    assert np.allclose(mol.h_vaporization, 40.66, atol=1e-2)
    assert np.allclose(mol.h_fusion, 13.3, atol=1e-1)
    assert np.allclose(mol.viscosity_liquid(333.8), 7.26e-4, atol=1e-6)
    assert np.allclose(mol.viscosity_liquid(374.4), 4.92e-4, atol=1e-6)
    assert np.allclose(mol.viscosity_liquid(403.1), 3.91e-4, atol=1e-6)
    assert np.allclose(mol.viscosity_liquid(423.3), 3.40e-4, atol=1e-6)


@pytest.mark.joback
def test_p_dichlorobenzene_real_nbt():
    mol = joback.get_groups(
        "C1=CC(=CC=C1Cl)Cl", "smiles", normal_boiling_point=447
    )
    assert np.allclose(mol.critical_temperature, 681, atol=1)


@pytest.mark.joback
def test_acentric_factor():
    # Perfluorohexane
    mol = joback.get_groups(
        "C(C(C(C(F)(F)F)(F)F)(F)F)(C(C(F)(F)F)(F)F)(F)F",
        "smiles",
        normal_boiling_point=329.8,
    )
    assert np.allclose(mol.acentric_factor, 0.525, atol=1e-2)


@pytest.mark.joback
def test_vapor_pressure():
    mol = joback.get_groups(
        "CC(C)O", identifier_type="smiles", normal_boiling_point=355.4
    )
    assert np.allclose(mol.vapor_pressure(450), 15.09, atol=2)


@pytest.mark.joback
def test_liquid_heat_capacity():
    mol = joback.get_groups("CC(=O)C", identifier_type="smiles")
    assert np.allclose(28.0, mol.heat_capacity_liquid(180) * 0.239, atol=4e-1)
    assert np.allclose(28.2, mol.heat_capacity_liquid(209) * 0.239, atol=1)
    assert np.allclose(29.8, mol.heat_capacity_liquid(297) * 0.239, atol=3)


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
    assert mol1.h_formation == mol2.h_formation
    assert mol1.g_formation == mol2.g_formation
    assert np.allclose(
        mol1.heat_capacity_ideal_gas_params,
        mol2.heat_capacity_ideal_gas_params,
    )
    assert mol1.h_fusion == mol2.h_fusion
    assert mol1.h_vaporization == mol2.h_vaporization
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
