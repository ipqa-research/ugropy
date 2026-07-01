import pytest

from ugropy import DefaultSolver


def test_ilp_solver():
    problem = DefaultSolver([], {})

    with pytest.raises(NotImplementedError):
        super(DefaultSolver, problem).solve_one_problem([])

    with pytest.raises(NotImplementedError):
        super(DefaultSolver, problem).solve()
