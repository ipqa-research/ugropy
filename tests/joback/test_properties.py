import numpy as np

from ugropy import Joback


def test_p_dichlorobenzene():
    mol = Joback("C1=CC(=CC=C1Cl)Cl", "smiles")

    assert mol.groups == {"-Cl": 2, "ring=CH-": 4, "ring=C<": 2}
    assert np.allclose(mol.normal_boiling_temperature, 443.4, atol=1e-2)
    assert np.allclose(mol.fusion_temperature, 256, atol=1)
    assert np.allclose(mol.critical_temperature, 675, atol=1)
    assert np.allclose(mol.critical_pressure, 41.5, atol=1e-1)
    assert np.allclose(mol.critical_volume, 362, atol=1)
    assert np.allclose(mol.h_formation, 26.41, atol=1e-2)
    assert np.allclose(mol.g_formation, 78.56, atol=1e-2)
    assert np.allclose(mol.heat_capacity(298), 112.3, atol=1)
    assert np.allclose(mol.heat_capacity(400), 139.2, atol=1)
    assert np.allclose(mol.heat_capacity(800), 206.8, atol=1)
    assert np.allclose(mol.heat_capacity(1000), 224.6, atol=1)
    assert np.allclose(mol.h_vaporization, 40.66, atol=1e-2)
    assert np.allclose(mol.h_fusion, 13.3, atol=1e-1)
    assert np.allclose(mol.liquid_viscosity(333.8), 7.26e-4, atol=1e-6)
    assert np.allclose(mol.liquid_viscosity(374.4), 4.92e-4, atol=1e-6)
    assert np.allclose(mol.liquid_viscosity(403.1), 3.91e-4, atol=1e-6)
    assert np.allclose(mol.liquid_viscosity(423.3), 3.40e-4, atol=1e-6)
