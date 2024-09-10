from typing import List

import itertools

from ugropy.refactor.fragment import Fragment

from rdkit import Chem

import numpy as np

import pulp


class FragmentationModel:
    def __init__(self, fragments: List[Fragment]):
        self.fragments = fragments

    def detect_fragments(self, molecule: Chem.rdchem.Mol):
        batch = DetectionBatch(molecule)

        for fragment in self.fragments:
            match = molecule.GetSubstructMatches(fragment.mol_object)

            if match:
                batch.add_fragment(fragment.name, match)
        return batch


class DetectionBatch:
    def __init__(self, molecule: Chem.rdchem.Mol):
        self.n = molecule.GetNumAtoms()
        self.fragments = {}
        self.overlaped_fragments = {}
        self.selected_fragments = []
        self.overlaped_atoms = []
        self.atoms_matrix = np.array([])
        self.solution_atoms = {}
        self.solution = {}

        self.has_overlap = False

    def get_groups(self):
        self.build_overlap_matrix()
        self.get_overlaped_fragments()
        self.solve_overlap()

        for frag in self.overlaped_fragments.keys():
            if frag not in self.selected_fragments:
                self.fragments.pop(frag)

        for frag in self.fragments.keys():
            name = frag.split("_")[0]

            if name not in self.solution_atoms.keys():
                self.solution_atoms[name] = [self.fragments[frag]]
            else:
                self.solution_atoms[name].append(self.fragments[frag])

        for frag in self.solution_atoms.keys():
            self.solution[frag] = len(self.solution_atoms[frag])

    def add_fragment(self, fragment_name: str, fragments: tuple):
        for i, f in enumerate(fragments):
            self.fragments[f"{fragment_name}_{i}"] = list(f)

    def build_overlap_matrix(self):
        self.atoms_matrix = np.zeros((len(self.fragments), self.n))

        for i, fragment in enumerate(self.fragments.values()):
            self.atoms_matrix[i, fragment] = 1

    def get_overlaped_fragments(self):
        overlap = np.sum(self.atoms_matrix, axis=0)
        self.overlaped_atoms = np.argwhere(overlap > 1).flatten()

        for name, frag in self.fragments.items():
            if np.isin(frag, self.overlaped_atoms).any():
                self.overlaped_fragments[name] = frag

    def solve_overlap(self):
        universe = set(self.overlaped_atoms)

        all_elements = set(
            itertools.chain.from_iterable(self.overlaped_fragments.values())
        )

        universe.update(all_elements)

        problem = pulp.LpProblem("Set_Cover_Problem", pulp.LpMinimize)

        n_frag = len(self.overlaped_fragments)

        x = pulp.LpVariable.dicts("x", range(n_frag), cat="Binary")

        problem += pulp.lpSum([x[i] for i in range(n_frag)])

        for elem in universe:
            sum_list = []
            for i, subset in enumerate(self.overlaped_fragments.values()):
                if elem in subset:
                    sum_list.append(x[i])

            # print(f"Restricción para el elemento {elem}: {sum_list} == 1")
            problem += pulp.lpSum(sum_list) == 1

        solver = pulp.getSolver("PULP_CBC_CMD", msg=False)

        problem.solve(solver)

        selected_subsets = [
            name
            for i, name in enumerate(self.overlaped_fragments.keys())
            if pulp.value(x[i]) == 1
        ]

        self.selected_fragments = selected_subsets
