"""ILPSolver abstract class.

Template for defining ILP solvers used by `ugropy` to solve the `Set Cover`
problem.
"""

import itertools
from abc import ABC, abstractmethod
from typing import List


class ILPSolver(ABC):
    """ILPSolver abstract class.

    Template for defining ILP solvers used by `ugropy` to solve the `Set Cover`
    problem.

    To inherit from this class, the following methods must be implemented:

    - `solve_one_problem`:
        Method that solves one instance of the `Set Cover` problem. This method
        must return a list of integers representing the chosen fragments and
        also must append the solution to the `selected_fragments` attribute
        as a dictionary.
    - `solve`:
        This class call the `solve_one_problem` method multiples times if the
        `search_multiple_solutions` attribute is True. The successives calls
        must include the previous solutions obtained

    You can personalize both methods as you want. The library will call the
    `solve` method and at the end, the `selected_fragments` attribute must
    contain the selected fragments as a list of dictionaries.

    Parameters
    ----------
    overlapped_atoms : List[int]
        List of atoms that are overlapped. On the `Set Cover` problem therms,
        these atoms represents the universe of the problem.
    fragments : dict[str: List[int]]
        Dictionary with the fragments (groups) to be selected. On the `Set
        Cover` problem therms, these fragments represents the sets to be
        selected to cover all the universe. The key is the name of the fragment
        and the value is a list of atoms that compose the fragment. Some atoms
        that represents a fragment could not be present in the universe. This
        is handled by adding thos "free atoms" to the universe, forcing that
        fragment to appear on the solution.
    search_multiple_solutions : bool, optional
        If True, the solver will search for multiple solutions, by default
        False

    Attributes
    ----------
    overlapped_atoms : List[int]
        List of atoms that are overlapped. On the `Set Cover` problem therms,
        these atoms represents the universe of the problem.
    fragments : dict[str: List[int]]
        Fragments to choose.
    search_multiple_solutions : bool
        If True, the solver will search for multiple solutions.
    selected_fragments : List[dict]
        Solutions found by the `solve` method.
    """

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
        """Solves one instance of the `Set Cover` problem.

        Parameters
        ----------
        not_valid_solutions : List[list[int]]
            List with binary values representing the fragments that are choosen
            with a 1 and the not choosen with a 0. This list contains the lists
            of the solutions already found.

        Returns
        -------
        List[int]
            Solution found by the solver. The list must contain binary values
            representing the fragments that are choosen with a 1 and the not
            choosen with a 0.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented
        """
        raise NotImplementedError(
            "Abstract method `solve_one_problem` not implemented"
        )

    @abstractmethod
    def solve(self) -> None:
        """Solve method called by fragmentations models.

        Ath the end of this method's execution, the `selected_fragments` must
        contain the selected fragments (solutions) as a list of dictionaries.

        Raises
        ------
        NotImplementedError
            Abstract method not implemented
        """
        raise NotImplementedError("Abstract method `solve` not implemented")
