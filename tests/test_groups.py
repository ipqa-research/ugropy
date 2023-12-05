import numpy as np

from rdkit import Chem

from ugropy import (
    Groups,
    get_joback_groups,
    get_psrk_groups,
    get_unifac_groups,
)


def test_smiles():
    mol = Groups("CCO", identifier_type="smiles")

    assert np.allclose(mol.molecular_weight, 46.07, atol=1e-2)

    assert mol.unifac_groups == get_unifac_groups(
        "CCO", identifier_type="smiles"
    )
    assert mol.psrk_groups == get_psrk_groups("CCO", identifier_type="smiles")
    assert mol.joback.groups == get_joback_groups(
        "CCO", identifier_type="smiles"
    )


def test_name():
    mol = Groups("ethanol", identifier_type="name")

    assert mol.unifac_groups == get_unifac_groups(
        "ethanol", identifier_type="name"
    )
    assert mol.psrk_groups == get_psrk_groups(
        "ethanol", identifier_type="name"
    )
    assert mol.joback.groups == get_joback_groups(
        "ethanol", identifier_type="name"
    )


def test_mol():
    chm = Chem.MolFromInchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    mol = Groups(chm, identifier_type="mol")

    assert mol.unifac_groups == get_unifac_groups(chm, identifier_type="mol")
    assert mol.psrk_groups == get_psrk_groups(chm, identifier_type="mol")
    assert mol.joback.groups == get_joback_groups(chm, identifier_type="mol")
