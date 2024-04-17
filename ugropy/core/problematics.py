"""Correct problematic structures module.

The algorithm of the function detect_groups may have some troubles with some
chemical structures. For example, in an ester group, it will detect also an
ether group. This problem can be handled by changing the SMARTS representation
of the ether group to something like:

[CH3][O][[#6]&!$([C]=O)]

With this SMARTS representation of the CH3O UNIFAC group, we are specifying
that this functional group it's a methyl group (CH3) bounded to oxygen by a
simple covalent bond, and that oxygen is bonded to any carbon but a carbon that
it's double bonded to an oxygen. This should avoid detecting an ether group in
an ester group. But, consider the structure of the molecule ethyl methyl
carbonate (PubChem CID 522046). That molecule has both an ester group and an
ether group, and the previous smarts representation will not detect the ether
group that we want to be detected. This problem defines what a problematic
structure is for the ugropy library. Maybe there is a SMARTS representation
that well behaves in these situations, but it's easier to make a list of these
problematic structures and correct them with a function than generate a complex
SMARTS representation of a functional group.
"""

import json

from rdkit import Chem

from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


def correct_problematics(
    mol_object: Chem.rdchem.Mol,
    mol_subgroups: dict,
    model: FragmentationModel,
) -> dict:
    """Correct problematic structures in mol_object.

    Parameters
    ----------
    mol_object : Chem.rdchem.Mol
        RDKit Chem object
    mol_subgroups : dict
        Dictionary with the subgroups not problematics corrected in mol_object.
    model: FragmentationModel
        FragmentationModel object.

    Returns
    -------
    dict
        Molecule's subrgoups corrected with the problematic structures list.
    """
    corrected_subgroups = mol_subgroups.copy()

    for smarts in model.problematic_structures.index:
        structure = Chem.MolFromSmarts(smarts)
        matches = mol_object.GetSubstructMatches(structure)
        how_many_problems = len(matches)

        if how_many_problems > 0:
            problematic_dict = json.loads(
                model.problematic_structures.loc[smarts, "contribute"]
            )

            for subgroup, contribution in problematic_dict.items():
                corrected_subgroups[subgroup] = (
                    corrected_subgroups.get(subgroup, 0)
                    + contribution * how_many_problems
                )

    # Eliminate occurrences == 0
    corrected_subgroups = {
        key: value for key, value in corrected_subgroups.items() if value != 0
    }
    return corrected_subgroups
