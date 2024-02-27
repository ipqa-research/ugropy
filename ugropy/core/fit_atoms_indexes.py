from itertools import combinations, product, chain

import numpy as np

from rdkit import Chem

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


def fit_atoms(
    mol_object: Chem.rdchem.Mol, mol_subgroups: dict, model: FragmentationModel
):
    # =========================================================================
    # Number of atoms in mol_object
    # =========================================================================
    total_atom_num = mol_object.GetNumAtoms()
    subgroups = list(mol_subgroups.keys())

    # =========================================================================
    # Getting true smarts, this will not be necessary later
    # =========================================================================
    smarts = {}
    for group in mol_subgroups.keys():
        s = model.subgroups.loc[group, "true_smarts"]

        if s is np.nan:
            s = model.subgroups.loc[group, "smarts"]

        smarts.update({group: s})

    # =========================================================================
    # Getting atoms candidates for each group
    # =========================================================================
    groups_atoms = {}
    for group, smt in smarts.items():
        atom = mol_object.GetSubstructMatches(Chem.MolFromSmarts(smt))
        groups_atoms[group] = atom

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
    # TODO: can be optimized? Probably no, but with an algorithm change yes
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
