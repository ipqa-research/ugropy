from typing import List

from ugropy.core.frag_classes.gibbs_model.gibbs_model import GibbsModel
from ugropy.core.frag_classes.gibbs_model.gibbs_result import (
    GibbsFragmentationResult,
)


def filter_big_solutions(
    solutions: List[GibbsFragmentationResult],
    gibbs_model: GibbsModel,
) -> List[GibbsFragmentationResult]:
    """Filter solutions based on their R value.

    Parameters
    ----------
    solutions : List[GibbsFragmentationResult]
        List of Gibbs fragmentation results to filter.
    gibbs_model : GibbsModel
        Gibbs model used to calculate molecular weights.

    Returns
    -------
    List[GibbsFragmentationResult]
        Filtered list of Gibbs fragmentation results.
    """
    if len(solutions) <= 1:
        return solutions

    # Calculate molecular weights for each solution
    mw_solutions = []
    for sol in solutions:
        mw = gibbs_model.calculate_molecular_weight(sol.subgroups)
        mw_solutions.append((mw, sol))

    # Find the minimum molecular weight
    min_mw = min(mw for mw, _ in mw_solutions)

    # Filter solutions that have the minimum molecular weight
    filtered_solutions = [sol for mw, sol in mw_solutions if mw == min_mw]

    return filtered_solutions
