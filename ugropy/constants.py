from pathlib import Path

import pandas as pd


# constants.py path
here = Path(__file__).parent.resolve()

# Dataframes

# UNIFAC
with open(f"{here}/groupscsv/unifac/unifac_subgroups.csv", mode='r') as f:
    unifac_subgroups = pd.read_csv(f, sep='|', index_col='group')

with open(f"{here}/groupscsv/unifac/unifac_matrix.csv", mode='r') as f:
    unifac_matrix = pd.read_csv(f, sep='|', index_col='group')

with open(f"{here}/groupscsv/unifac/ch2_hideouts.csv", mode='r') as f:
    ch2_hideouts = pd.read_csv(f, index_col='group').index
    
with open(f"{here}/groupscsv/unifac/ch_hideouts.csv", mode='r') as f:
    ch_hideouts = pd.read_csv(f, index_col='group').index

# Problematics
with open(f"{here}/groupscsv/problematic_structures.csv", mode='r') as f:
    problematic_structures = pd.read_csv(f, sep='|', index_col='smarts', comment="?")
