"""fit_atoms_indexes module."""

from itertools import chain, combinations, product

from rdkit import Chem

from ugropy.core.detect_model_groups import group_matches
from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


def fit_atoms(
    mol_object: Chem.rdchem.Mol, mol_subgroups: dict, model: FragmentationModel
) -> dict:
    """Assign the atoms indexes for each mol_subgroup.

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Mol object.
    mol_subgroups : dict
        Subgroups of mol_object.
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    dict
        Atom indexes in mol_object of each subgroup.
    """
    # =========================================================================
    # Number of atoms in mol_object and subgroups.
    # =========================================================================
    total_atom_num = mol_object.GetNumAtoms()
    subgroups = list(mol_subgroups.keys())

    # import ipdb; ipdb.set_trace()

    # =========================================================================
    # Getting atoms candidates for each group.
    # =========================================================================
    groups_atoms = {}
    for group in mol_subgroups.keys():
        atoms = group_matches(mol_object, group, model, "fit")
        groups_atoms[group] = atoms

    # =========================================================================
    # Getting combinations for each subgroup according to the number of the
    # number of occurences in the mol_subgroups tentative solution.
    # =========================================================================
    atoms_combinations = {
        group: list(combinations(groups_atoms[group], mol_subgroups[group]))
        for group in mol_subgroups
    }

    # =========================================================================
    # Check all possible combinations of solutions
    # TODO: can be optimized? Probably not, but with an algorithm change maybe
    # =========================================================================
    for comb in product(*atoms_combinations.values()):
        plain_tuple = chain(*chain(*comb))

        if len(set(plain_tuple)) == total_atom_num:
            sol_comb_dict = {
                subgroups[i]: (comb[i][0] if len(comb[i]) == 1 else comb[i])
                for i in range(len(subgroups))
            }
            return sol_comb_dict

    return {}
