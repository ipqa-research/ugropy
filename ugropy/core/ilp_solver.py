from abc import ABC, abstractmethod
from typing import List


class ILPSolver(ABC):
    def __init__(self, overlapped_atoms: List, overlapped_fragments: dict):
        self.overlapped_atoms = overlapped_atoms
        self.overlapped_fragments = overlapped_fragments
        self.selected_fragments = {}

    @abstractmethod
    def solve(self):
        raise NotImplementedError
