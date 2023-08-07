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
    df_problematics = ...

    # Groups and occurence in name:
    groups = np.array([])
    many_groups = np.array([])

    for group in df.index:
        try:
            #import ipdb
            #ipdb.set_trace(cond=group == "CH2O")
            
            group_count = 0
            smarts = df.loc[group]["smarts"]

            func_group = Chem.MolFromSmarts(smarts)
            matches = chem_object.GetSubstructMatches(func_group)
            how_many = len(matches)
            
            if how_many > 0:
                group_count += how_many
                if group not in groups:
                    groups = np.append(groups, group)
            
            if group_count > 0:
                many_groups = np.append(many_groups, group_count).astype(int)
                #import ipdb
                #ipdb.set_trace(cond= (many_groups == np.array([3, 2, 1, 3, 2])).all())

        except:
            ...

    #import ipdb
    #ipdb.set_trace()

    # get final dict
    dff = dfm.loc[groups][groups]
    dff = dff.mul(many_groups, axis= 0)
    dff_sum = dff.sum(axis=0)
    dff_sum.replace(0, pd.NA, inplace=True)
    dff_final = dff_sum.dropna()

    return dff_final.to_dict()
