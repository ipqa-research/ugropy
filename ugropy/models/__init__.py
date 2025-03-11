"""Fragmentation models implementations.

All models can be imported directly with:

.. code-block:: python

    from ugropy import joback, psrk, unifac, abdulelah_gani

You can check the group list and their SMARTS representation with:

.. code-block:: python

    joback.subgroups
    psrk.subgroups
    unifac.subgroups

Obtain the contribution of each group to the properties with:

.. code-block:: python

    joback.get_groups("n-hexane").subgroups
    unifac.get_groups("toluene").subgroups
    psrk.get_groups("CCO", "smiles").subgroups
    abdulelah_gani.get_groups("CCO", "smiles").primary.subgroups
    abdulelah_gani.get_groups("CCO", "smiles").secondary.subgroups
    abdulelah_gani.get_groups("CCO", "smiles").tertiary.subgroups

Attributes
----------
psrk: GibbsModel
    Predictive Soave-Redlich-Kwong FragmentationModel :cite:p:`ddbst, psrk1,
    psrk2`
unifac : GibbsModel
    Classic LV-UNIFAC FragmentationModel :cite:p:`ddbst, unifac1, unifac2,
    unifac3, unifac4, unifac5, unifac6`
dortmund : GibbsModel
    Dortmund mod-UNIFAC FragmentationModel :cite:p:`dor1, dor2, dor3, dor4,
    dor5, dor6, dor7, dor8, dor9, dor10, dor11`
joback: JobackModel
    Joback FragmentationModel :cite:p:`joback1, joback2`
abdulelah_gani_p: AbdulelahGaniPSTModel
    Abdulelah-Gani Primary Structures FragmentationModel :cite:p:`gani`
abdulelah_gani_s: AbdulelahGaniPSTModel
    Abdulelah-Gani Secondary Structures FragmentationModel :cite:p:`gani`
abdulelah_gani_t: AbdulelahGaniPSTModel
    Abdulelah-Gani Tertiary Structures FragmentationModel :cite:p:`gani`
abdulelah_gani: AbdulelahGaniModel
    Abdulelah-Gani properties estimation model :cite:p:`gani`
"""

from . import (
    abdulelah_gani_pmod,
    abdulelah_gani_smod,
    abdulelah_gani_tmod,
    dortmundmod,
    jobackmod,
    psrkmod,
    read_csv,
    unifacmod,
)


__all__ = [
    "abdulelah_gani_pmod",
    "abdulelah_gani_smod",
    "abdulelah_gani_tmod",
    "dortmundmod",
    "read_csv",
    "jobackmod",
    "psrkmod",
    "unifacmod",
]
