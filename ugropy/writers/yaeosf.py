"""to_yaeos module."""

from typing import List

from ugropy.core.frag_classes.gibbs_model.gibbs_model import GibbsModel
from ugropy.writers.thermo import to_thermo


def to_yaeos(mol_subgroups_list: List[dict], model: GibbsModel) -> str:
    """Obtain the Fortran source code for yaeos groups definition.

    yaeos: https://github.com/ipqa-research/yaeos

    Parameters
    ----------
    mol_subgroups_list : List[dict]
        List of ugropy subgroups dictionaries.
    model : GibbsModel
        Gibbs excess FragmentationModel (unifac, psrk, etc).

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

    # Obtain the dictionary of subgroups for each molecule but with index
    for i, mol_subgroups in enumerate(mol_subgroups_list):
        thermo_groups = to_thermo(mol_subgroups, model)

        molecule_code = (
            f"molecules({i+1})%groups_ids = "
            f"{list(thermo_groups.keys())}\n"
            f"molecules({i+1})%number_of_groups = "
            f"{list(thermo_groups.values())}\n"
            "\n"
        )

        code += molecule_code

    return code
