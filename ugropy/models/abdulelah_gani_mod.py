from ugropy.constants import _csvs
from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani import AbdulelahGaniModel
from ugropy.models.abdulelah_gani_pmod import abdulelah_gani_p
from ugropy.models.abdulelah_gani_smod import abdulelah_gani_s
from ugropy.models.abdulelah_gani_tmod import abdulelah_gani_t


abdulelah_gani = AbdulelahGaniModel(abdulelah_gani_p, abdulelah_gani_s, abdulelah_gani_t)
