"""ilp_solvers module.

This module contains the integer linear programming solvers used by `ugropy`.
ILP solvers are neccesary to solve the `Set Cover` problem, which is solved
to select the optimal set of fragments (groups) to cover the overlapping atoms
on a molecule.
"""

from . import default_solver, ilp_solver


__all__ = ["ilp_solver", "default_solver"]
