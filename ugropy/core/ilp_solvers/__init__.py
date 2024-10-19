"""ilp_solvers module.

This module contains the integer linear programming solvers used by `ugropy`.
ILP solvers are neccesary to solve the `Set Cover` problem, which is solved
to select the optimal set of fragments (groups) to cover the overlapping atoms
on a molecule.
"""

from .default_solver import DefaultSolver
from .ilp_solver import ILPSolver


__all__ = ["ILPSolver", "DefaultSolver"]
