"""to_thermo module."""

from ugropy.core.frag_classes.gibbs_model.gibbs_model import GibbsModel


def to_thermo(mol_subgroups: dict, model: GibbsModel) -> dict:
    """Obtain the subgroups dictionary to the Caleb Bell's Thermo library.

    Thermo: https://github.com/CalebBell/thermo

    Parameters
    ----------
    mol_subgroups : dict
        ugropy subgroups.
    model : GibbsModel
        Gibbs excess FragmentationModel (unifac or psrk).

    Returns
    -------
    dict
        Thermo fragmentation subgroups.
    """
    thermo_groups = {}
    for group, occurrence in mol_subgroups.items():
        group_num = model.subgroups_info.loc[group, "subgroup_number"]

        thermo_groups[group_num] = occurrence

    return thermo_groups
