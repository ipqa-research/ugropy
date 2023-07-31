from rdkit import Chem

import pubchempy as pcp

import numpy as np

import pandas as pd


def get_groups(name, subgroups, subgroups_matrix):
    pcp_object = pcp.get_compounds(name, "name")[0]
    smiles = pcp_object.canonical_smiles
    chem_object = Chem.MolFromSmiles(smiles)

    df = subgroups
    dfm = subgroups_matrix

    # Groups and occurence in name:
    groups = np.array([])
    many_groups = np.array([])

    for group in df.index:
        try:
            smarts = df.loc[group]["smarts"].split()

            for idx, s in enumerate(smarts):
                count = 0

                func_group = Chem.MolFromSmarts(s)
                matches = chem_object.GetSubstructMatches(func_group)
                how_many = len(matches)
                
                if how_many > 0:
                    count += how_many
                    if idx == 0:
                        groups = np.append(groups, group)
            
            if count > 0:
                many_groups = np.append(many_groups, count).astype(int)
        except:
            ...

    # get final dict
    dff = dfm.loc[groups][groups]
    dff = dff * many_groups
    dff_sum = dff.sum(axis=0)
    dff_sum.replace(0, pd.NA, inplace=True)
    final_dict = dff_sum.to_dict()

    return final_dict
