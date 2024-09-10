from rdkit import Chem


class Fragment:
    def __init__(self, name: str, smarts: str):
        self.name = name
        self.smarts = smarts
        self.mol_object = Chem.MolFromSmarts(smarts)
