from ugropy.constants import _csvs
from ugropy.core.frag_classes.gibss_model.gibbs_model import GibbsModel
from ugropy.models.read_csv import _rd


# =============================================================================
# PSRK
# =============================================================================
_psrk = _csvs / "psrk"

_psrk_sg = _rd(_psrk / "psrk_subgroups.csv", "group")
_psrk_info = _rd(_psrk / "psrk_info.csv", "group")

psrk = GibbsModel(subgroups=_psrk_sg, subgroups_info=_psrk_info)
