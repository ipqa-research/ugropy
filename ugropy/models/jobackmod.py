"""Joback FragmentationModel implementation.

Import and use the Joback FragmentationModel with:

.. code-block:: python

    from ugropy import joback

    # Get groups from molecule's name
    tol = joback.get_groups("toluene")

    print(tol.subgroups)

    # Get groups from molecule's SMILES
    eth = joback.get_groups("CCO", "smiles")

    print(eth.subgroups)

Attributes
----------
joback: PropertiesEstimator
    Joback FragmentationModel :cite:p:`joback1, joback2`
"""

from ugropy.constants import _csvs
from ugropy.core.frag_classes.joback.joback_model import JobackModel
from ugropy.models.read_csv import _rd

# =============================================================================
# Joback
# =============================================================================
_jb = _csvs / "joback"

_jb_sg = _rd(_jb / "joback_subgroups.csv", "group")
_jb_props = _rd(_jb / "properties_contrib.csv", "group")

joback = JobackModel(_jb_sg, _jb_props)
