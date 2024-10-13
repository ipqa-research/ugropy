from typing import List

import pulp

import numpy as np

from ugropy.core.ilp_solvers.ilp_solver import ILPSolver


class DefaultSolver(ILPSolver):
    def solve_one_problem(self, not_valid_solutions: List = []) -> List[int]:
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

        # Constaints 2: if are not_valid_solutions, the newones cannot use less
        # elements than the previous ones, otherwise the first one was not
        # optimal.
        if not_valid_solutions:
            opt = np.sum(np.array(not_valid_solutions[0], dtype=int))

            problem += pulp.lpSum([x[i] for i in range(n_frag)]) == opt

        # Constaints 3: not already get solution
        for solution in not_valid_solutions:
            problem += (
                pulp.lpSum([x[i] * solution[i] for i in range(n_frag)])
                <= opt - 1
            )

        # Solver configuration verbosity
        solver = pulp.getSolver("PULP_CBC_CMD", msg=False)

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
            return None

    def solve(self) -> None:
        # =====================================================================
        # Multiple solutions
        # =====================================================================
        if self.search_multiple_solutions:
            not_valid_solutions = []
            while True:
                solution = self.solve_one_problem(not_valid_solutions)
                if solution is None:
                    break
                not_valid_solutions.append(solution)
        else:
            self.solve_one_problem()
        return self.selected_fragments
