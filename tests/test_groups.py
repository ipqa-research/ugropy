import numpy as np

from rdkit import Chem

from ugropy import Groups, joback, psrk, unifac


def test_smiles():
    mol = Groups("CCO", identifier_type="smiles")

    assert np.allclose(mol.molecular_weight, 46.07, atol=1e-2)

    assert (
        mol.unifac.subgroups
        == unifac.get_groups("CCO", identifier_type="smiles").subgroups
    )
    assert (
        mol.psrk.subgroups
        == psrk.get_groups("CCO", identifier_type="smiles").subgroups
    )
    assert (
        mol.joback.subgroups
        == joback.get_groups("CCO", identifier_type="smiles").subgroups
    )


def test_name():
    mol = Groups("ethanol", identifier_type="name")

    assert (
        mol.unifac.subgroups
        == unifac.get_groups("ethanol", identifier_type="name").subgroups
    )
    assert (
        mol.psrk.subgroups
        == psrk.get_groups("ethanol", identifier_type="name").subgroups
    )
    assert (
        mol.joback.subgroups
        == joback.get_groups("ethanol", identifier_type="name").subgroups
    )


def test_mol():
    chm = Chem.MolFromInchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    mol = Groups(chm, identifier_type="mol")

    assert (
        mol.unifac.subgroups
        == unifac.get_groups(chm, identifier_type="mol").subgroups
    )
    assert (
        mol.psrk.subgroups
        == psrk.get_groups(chm, identifier_type="mol").subgroups
    )
    assert (
        mol.joback.subgroups
        == joback.get_groups(chm, identifier_type="mol").subgroups
    )
