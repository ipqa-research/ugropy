from ugropy.refactor.fragment import Fragment
from ugropy.refactor.fragmentation_model import FragmentationModel


_unifac_fragments = [
    Fragment("CH3", "[CX4H3]"),
    Fragment("CH2", "[CX4H2]"),
    Fragment("CH", "[CX4H]"),
    Fragment("C", "[CX4H0]"),
    Fragment("OH", "[OH]"),
    Fragment("H2O", "[OH2]"),
    Fragment("ACH", "[cH]"),
    Fragment("AC", "[cH0]"),
    Fragment("ACCH3", "[cH0][CX4H3]"),
    Fragment("ACCH2", "[cH0][CX4H2]"),
    Fragment("ACCH", "[cH0][CX4H]"),
    Fragment("ACOH", "[cH0][OH]"),
    Fragment("CH3O", "[CH3][OH0]"),
    Fragment("CH2O", "[CH2][OH0]"),
    Fragment("CHO", "[CH][OH0]"),
    Fragment("COO", "[CX3H0](=[OH0])[OH0]"),
    Fragment("ACNH2", "[cH0][NH2]"),
]


unifac2 = FragmentationModel(_unifac_fragments)
