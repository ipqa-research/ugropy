from ugropy.constants import _csvs
from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani import (
    AbdulelahGaniModel,
)
from ugropy.models.abdulelah_gani_pmod import abdulelah_gani_p
from ugropy.models.abdulelah_gani_smod import abdulelah_gani_s
from ugropy.models.abdulelah_gani_tmod import abdulelah_gani_t
from ugropy.models.read_csv import _rd

# =============================================================================
# Abdulelah Gani Model
# =============================================================================
_ag = _csvs / "abdulelah_gani"

_ag_pc = _rd(_ag / "properties_contributions.csv", "group")
_ag_biases = _rd(_ag / "properties_biases.csv", "biases")


abdulelah_gani = AbdulelahGaniModel(
    abdulelah_gani_p,
    abdulelah_gani_s,
    abdulelah_gani_t,
    _ag_pc,
    _ag_biases,
)
