"""to_yaeos module."""

from typing import List

from ugropy.core.frag_classes.gibbs_model.gibbs_result import (
    GibbsFragmentationResult,
)


def to_yaeos(mol_subgroups_list: List[GibbsFragmentationResult]) -> str:
    """Obtain the Fortran source code for yaeos groups definition.

    yaeos: https://github.com/ipqa-research/yaeos

    Parameters
    ----------
    mol_subgroups_list : List[GibbsFragmentationResult]
        List of ugropy GibbsModel solutions (UNIFAC, PSRK, Dortmund, etc).

    Returns
    -------
    str
        Yaeos Fortran source code.
    """
    n_mol = len(mol_subgroups_list)

    code = (
        "use yaeos__models_ge_group_contribution_unifac, only: Groups\n"
        "\n"
        f"type(Groups) :: molecules({n_mol})\n"
        "\n"
    )

    subgroups_numbers = [g.subgroups_num for g in mol_subgroups_list]

    for i, mol_subgroups in enumerate(subgroups_numbers):
        molecule_code = (
            f"molecules({i+1})%groups_ids = "
            f"{list(mol_subgroups.keys())}\n"
            f"molecules({i+1})%number_of_groups = "
            f"{list(mol_subgroups.values())}\n"
            "\n"
        )

        code += molecule_code

    return code
