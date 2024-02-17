import numpy as np

from rdkit import Chem

from ugropy import Groups, get_groups, joback, psrk, unifac


def test_smiles():
    mol = Groups("CCO", identifier_type="smiles")

    assert np.allclose(mol.molecular_weight, 46.07, atol=1e-2)

    assert mol.unifac_groups == get_groups(
        unifac, "CCO", identifier_type="smiles"
    )
    assert mol.psrk_groups == get_groups(psrk, "CCO", identifier_type="smiles")
    assert mol.joback.groups == get_groups(
        joback, "CCO", identifier_type="smiles"
    )


def test_name():
    mol = Groups("ethanol", identifier_type="name")

    assert mol.unifac_groups == get_groups(
        unifac, "ethanol", identifier_type="name"
    )
    assert mol.psrk_groups == get_groups(
        psrk, "ethanol", identifier_type="name"
    )
    assert mol.joback.groups == get_groups(
        joback, "ethanol", identifier_type="name"
    )


def test_mol():
    chm = Chem.MolFromInchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
    mol = Groups(chm, identifier_type="mol")

    assert mol.unifac_groups == get_groups(unifac, chm, identifier_type="mol")
    assert mol.psrk_groups == get_groups(psrk, chm, identifier_type="mol")
    assert mol.joback.groups == get_groups(joback, chm, identifier_type="mol")
