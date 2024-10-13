import itertools

from abc import ABC, abstractmethod
from typing import List


class ILPSolver(ABC):
    def __init__(
        self,
        overlapped_atoms: List,
        fragments: dict,
        search_multiple_solutions: bool = False,
    ):
        self.overlapped_atoms = overlapped_atoms
        self.fragments = fragments
        self.search_multiple_solutions = search_multiple_solutions
        self.selected_fragments = []

        # Get overlapped fragments (which participates on overlapped atoms)
        self.overlapped_fragments = {}

        for name, fragment in self.fragments.items():
            if any(atom in self.overlapped_atoms for atom in fragment):
                self.overlapped_fragments[name] = fragment

        # set universe
        self.universe = set(self.overlapped_atoms)

        # Update universe with all elements in the overlapped fragments that
        # not present in the overlapped atoms but participates in a fragment.
        all_elements = set(
            itertools.chain.from_iterable(self.overlapped_fragments.values())
        )

        self.universe.update(all_elements)

    @abstractmethod
    def solve_one_problem(self, not_valid_solutions: List) -> List[int]:
        raise NotImplementedError(
            "Abstract method `solve_one_problem` not implemented"
        )

    @abstractmethod
    def solve(self) -> None:
        raise NotImplementedError("Abstract method `solve` not implemented")
