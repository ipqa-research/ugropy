import itertools

import pulp

from ugropy.core.ilp_solver import ILPSolver


class DefaultSolver(ILPSolver):
    def solve(self):
        # Universe of elements
        universe = set(self.overlapped_atoms)

        # Update universe with all elements in the overlapped fragments that
        # not present in the overlapped atoms but participates in a fragment.
        all_elements = set(
            itertools.chain.from_iterable(self.overlapped_fragments.values())
        )

        universe.update(all_elements)

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
        for elem in universe:
            sum_list = []
            for i, subset in enumerate(self.overlapped_fragments.values()):
                if elem in subset:
                    sum_list.append(x[i])

            problem += pulp.lpSum(sum_list) == 1

        # Solver configuration verbosity
        solver = pulp.getSolver("PULP_CBC_CMD", msg=False)

        # Solve
        problem.solve(solver)

        # Selected fragments (solution)
        selected_subsets = [
            name
            for i, name in enumerate(self.overlapped_fragments.keys())
            if pulp.value(x[i]) == 1
        ]

        self.selected_fragments = selected_subsets
