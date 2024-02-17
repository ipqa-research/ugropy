"""Fragmentation models implemented.

Attributes
----------
unifac : FragmentationModel
    Classic LV-UNIFAC FragmentationModel
psrk: FragmentationModel
    Predictive Soave-Redlich-Kwong FragmentationModel
dortmund: FragmentationModel
    Dortmund UNIFAC FragmentationModel
joback: FragmentationModel
    Joback FragmentationModel
"""

from ugropy.constants import (
    dort_ch2_hide,
    dort_ch_hide,
    dort_maingroups,
    dort_problem,
    dort_subgroups,
    joback_problem,
    joback_subgroups,
    psrk_ch2_hide,
    psrk_ch_hide,
    psrk_maingroups,
    psrk_problem,
    psrk_subgroups,
    unifac_ch2_hide,
    unifac_ch_hide,
    unifac_maingroups,
    unifac_problem,
    unifac_subgroups,
)
from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


# LV-UNIFAC
unifac = FragmentationModel(
    subgroups=unifac_subgroups,
    main_groups=unifac_maingroups,
    problematic_structures=unifac_problem,
    ch2_hideouts=unifac_ch2_hide.to_list(),
    ch_hideouts=unifac_ch_hide.to_list(),
)

# PSRK
psrk = FragmentationModel(
    subgroups=psrk_subgroups,
    main_groups=psrk_maingroups,
    problematic_structures=psrk_problem,
    ch2_hideouts=psrk_ch2_hide.to_list(),
    ch_hideouts=psrk_ch_hide.to_list(),
)

# Dortmund
dortmund = FragmentationModel(
    subgroups=dort_subgroups,
    main_groups=dort_maingroups,
    problematic_structures=dort_problem,
    ch2_hideouts=dort_ch2_hide.to_list(),
    ch_hideouts=dort_ch_hide.to_list(),
)

# Joback
joback = FragmentationModel(
    subgroups=joback_subgroups, problematic_structures=joback_problem
)
