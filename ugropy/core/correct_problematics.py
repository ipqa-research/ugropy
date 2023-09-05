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

import pandas as pd

from rdkit import Chem


def correct_problematics(
    chem_object: Chem.rdchem.Mol,
    filtered_subgroups: pd.DataFrame,
    problematic_structures: pd.DataFrame,
) -> pd.DataFrame:
    """Correct problematic structures in chem_object.

    Parameters
    ----------
    chem_object : Chem.rdchem.Mol
        RDKit Chem object
    filtered_subgroups : pd.DataFrame
        Subgroups matrix DataFrame filtered for groups present in chem_object.
    problematic_structures : pd.DataFrame
        Problematic structures DataFrame.

    Returns
    -------
    pd.DataFrame
        Subgroups matrix DataFrame filtered with corrected problematic
        structures.
    """
    dff = filtered_subgroups.copy()

    for smarts in problematic_structures.index:
        structure = Chem.MolFromSmarts(smarts)
        matches = chem_object.GetSubstructMatches(structure)
        how_many_problems = len(matches)

        if how_many_problems > 0:
            p_dict = json.loads(problematic_structures.loc[smarts].contribute)
            for grp in p_dict.keys():
                try:
                    dff.loc[grp, grp] += p_dict[grp] * how_many_problems
                except KeyError:
                    dff[grp] = 0
                    dff.loc[grp] = 0
                    dff.loc[grp, grp] += p_dict[grp] * how_many_problems
    return dff
