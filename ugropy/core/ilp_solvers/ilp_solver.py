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
        must return a list of integers representing the chosen fragments. Also
        must append the solution to the `selected_fragments` attribute as a
        dictionary.
    - `solve`:
        This class call the `solve_one_problem` method multiples times if the
        `search_multiple_solutions` attribute is True. The successives calls
        must include the previous solutions obtained as restrictions.

    You can personalize both methods as you want. The library will call the
    `solve` method and at the end of its execution the `selected_fragments`
    attribute must contain the selected fragments for each solution.

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
        is handled by adding those "free atoms" to the universe, forcing that
        fragment to appear on the solution. This updating of the universe is
        done on the `__init__` method, so you don't need to worry about it.
    search_multiple_solutions : bool, optional
        If True, the solver will search for multiple solutions, by default
        False
    search_nonoptimal : bool, optional
        If True, the solver will search for non-optimal solutions along with
        the optimal ones. This is useful when the user wants to find all
        possible combinations of fragments that cover the universe. By default
        False. If `search_multiple_solutions` is True, this parameter will be
        ignored.
    solver_arguments : dict, optional
        Dictionary with the arguments to be passed to the solver.

    Attributes
    ----------
    overlapped_atoms : List[int]
        List of atoms that are overlapped. On the `Set Cover` problem therms,
        these atoms represents the universe of the problem.
    fragments : dict[str: List[int]]
        Fragments to choose to cover the universe.
    search_multiple_solutions : bool
        If True, the solver will search for multiple solutions.
    search_nonoptimal : bool
        If True, the solver will search for non-optimal solutions along with
        the optimal ones.
    selected_fragments : List[List[str]]
        Solutions found by the `solve` method. Its a list  that contains the
        different obtained solutions as lists of strings. Each string
        corresponds to the name of the selected fragments.
    overlapped_fragments : dict
        Fragments that participates on overlapped atoms.
    universe : set
        Set with all elements that must be covered by the fragments. This set
        is composed by the overlapped atoms and the atoms that are present in
        the overlapped fragments but not in the overlapped atoms.
    """

    def __init__(
        self,
        overlapped_atoms: List,
        fragments: dict,
        search_multiple_solutions: bool = False,
        search_nonoptimal: bool = False,
        solver_arguments: dict = {},
    ):
        self.overlapped_atoms = overlapped_atoms
        self.fragments = fragments
        self.search_multiple_solutions = search_multiple_solutions
        self.search_nonoptimal = search_nonoptimal

        if solver_arguments == {}:
            self.solver_arguments = {
                "solver": "PULP_CBC_CMD",
            }
        else:
            self.solver_arguments = solver_arguments

        # Attribute that will store the selected fragments after a
        # `solve` method call.
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
            choosen with a 0. Of course, the length of the list must be equal
            to the lenght of the `overlapped_fragments` attribute.

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
