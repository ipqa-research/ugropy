from typing import List

from ugropy.refactor.fragment import Fragment

from rdkit import Chem


class FragmentationModel:
    def __init__(self, fragments: List[Fragment]):
        self.fragments = fragments

    def detect_fragments(self, molecule: Chem.rdchem.Mol):
        detected = {}

        for fragment in self.fragments:
            match = molecule.GetSubstructMatches(fragment.mol_object)

            if match:
                detected[fragment.name] = match

        return detected