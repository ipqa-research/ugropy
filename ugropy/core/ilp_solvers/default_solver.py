"""Default ILP solver for the set cover problem.

The DefaultSolver uses the Pulp library to solve the set cover problem.
"""

from typing import List

import numpy as np

import pulp

from ugropy.core.ilp_solvers.ilp_solver import ILPSolver


class DefaultSolver(ILPSolver):
    """Default ILP solver for the set cover problem."""

    def solve_one_problem(self, not_valid_solutions: List = []) -> List[int]:
        """Solve one `Set Cover` problem with the Pulp library.

        Parameters
        ----------
        not_valid_solutions : List, optional
            Already got solutions, by default []

        Returns
        -------
        List[int]
            List of binary values representing the selected fragments.
        """
        # =====================================================================
        # ILP formulation (Pulp)
        # =====================================================================
        problem = pulp.LpProblem("Set_Cover_Problem", pulp.LpMinimize)

        n_frag = len(self.overlapped_fragments)

        # Decision variables
        x = pulp.LpVariable.dicts("x", range(n_frag), cat="Binary")

        # Objective function
        problem += pulp.lpSum([x[i] for i in range(n_frag)])

        # Constraints 1: each element must be covered by only one subset
        for elem in self.universe:
            sum_list = []
            for i, subset in enumerate(self.overlapped_fragments.values()):
                if elem in subset:
                    sum_list.append(x[i])

            problem += pulp.lpSum(sum_list) == 1

        # Constaints 2: if are not_valid_solutions, the new ones cannot use
        # less elements than the previous ones, otherwise the first one was not
        # optimal.
        if not_valid_solutions:
            opt = np.sum(np.array(not_valid_solutions[0], dtype=int))

            if not self.search_nonoptimal:
                # Do you want to search only optimals solutions?
                problem += pulp.lpSum([x[i] for i in range(n_frag)]) == opt

        # Constaints 3: already got solutions not allowed
        for solution in not_valid_solutions:
            problem += (
                pulp.lpSum([x[i] * solution[i] for i in range(n_frag)])
                <= pulp.lpSum([solution[i] for i in range(n_frag)]) - 1
            )

        # Solver configuration verbosity
        solver = pulp.getSolver(self.solver_arguments["solver"], msg=False)

        # Solve
        problem.solve(solver)

        # Selected fragments (solution)
        if pulp.LpStatus[problem.status] == "Optimal":
            selected_subsets = [
                name
                for i, name in enumerate(self.overlapped_fragments.keys())
                if pulp.value(x[i]) == 1
            ]

            self.selected_fragments.append(selected_subsets)

            return [pulp.value(x[i]) for i in range(n_frag)]
        else:
            # No feasible problem
            return None

    def solve(self) -> None:
        """Run multiple times the solve_one_problem method."""
        # Multiple solutions
        if self.search_multiple_solutions:
            not_valid_solutions = []
            while True:
                solution = self.solve_one_problem(not_valid_solutions)
                if solution is None:
                    break
                not_valid_solutions.append(solution)
        else:
            _ = self.solve_one_problem()
