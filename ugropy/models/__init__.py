"""Fragmentation models implementations.

All models can be imported directly with:

.. code-block:: python

    from ugropy import joback, psrk, unifac

You can check the group list and their SMARTS representation with:

.. code-block:: python

    joback.subgroups psrk.subgroups unifac.subgroups

Obtain the contribution of each group to the properties with:

.. code-block:: python

    joback.get_groups("n-hexane").subgroups
    unifac.get_groups("toluene").subgroups
    psrk.get_groups("CCO", "smiles").subgroups

Attributes
----------
psrk: GibbsModel
    Predictive Soave-Redlich-Kwong FragmentationModel :cite:p:`ddbst, psrk1,
    psrk2`
unifac : GibbsModel
    Classic LV-UNIFAC FragmentationModel :cite:p:`ddbst, unifac1, unifac2,
    unifac3, unifac4, unifac5, unifac6`
joback: JobackModel
    Joback FragmentationModel :cite:p:`joback1, joback2`
abdulelah_gani_p: AbdulelahGaniPrimaryModel
    Abdulelah-Gani Primary Structures FragmentationModel :cite:p:`gani`
"""

from . import abdulelah_gani_pmod, jobackmod, psrkmod, read_csv, unifacmod


__all__ = [
    "abdulelah_gani_pmod",
    "read_csv",
    "jobackmod",
    "psrkmod",
    "unifacmod",
]
