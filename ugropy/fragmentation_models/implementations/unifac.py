"""Fragmentation models implemented.

All models can be imported directly with:

.. code-block:: python

   from ugropy import joback, psrk, unifac

You can check the group list and their SMARTS representation with:

.. code-block:: python

   joback.subgroups psrk.subgroups unifac.subgroups

In the case of a PropertiesEstimator like joback, you can check the
contribution of each group to the properties with:

.. code-block:: python

    joback.properties_contributions

Attributes
----------
unifac : GibbsModel
    Classic LV-UNIFAC FragmentationModel :cite:p:`ddbst, unifac1, unifac2,
    unifac3, unifac4, unifac5, unifac6`
psrk: GibbsModel
    Predictive Soave-Redlich-Kwong FragmentationModel :cite:p:`ddbst, psrk1,
    psrk2`
joback: PropertiesEstimator
    Joback FragmentationModel :cite:p:`joback1, joback2`
"""

from ugropy.constants import _csvs
from ugropy.fragmentation_models.gibbs_model import GibbsModel
from ugropy.fragmentation_models.implementations.read_csv import _rd


# =============================================================================
# LV-UNIFAC
# =============================================================================
_uni = _csvs / "unifac"

_uni_sg = _rd(_uni / "unifac_subgroups.csv", "group")
_uni_info = _rd(_uni / "unifac_info.csv", "group")

unifac = GibbsModel(subgroups=_uni_sg, subgroups_info=_uni_info)
