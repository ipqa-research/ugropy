from ugropy.refactor.fragment import Fragment
from ugropy.refactor.fragmentation_model import FragmentationModel


_unifac_fragments = [
    Fragment("CH3", "[CH3]"),
    Fragment("CH2", "[CH2]"),
    Fragment("CH", "[CH]"),
    Fragment("C", "[CH0]"),
    Fragment("OH", "[OH]"),
    Fragment("H2O", "O"),
    Fragment("ACH", "[cH]"),
    Fragment("AC", "[cH0]"),
    Fragment("ACCH3", "[cH0][CH3]"),
    Fragment("ACCH2", "[cH0][CH2]"),
    Fragment("ACCH", "[cH0][CH]"),
]


unifac = FragmentationModel(_unifac_fragments)

# mol = Chem.MolFromSmiles("C(C1=CC=CC=C1)C1=CC=CC=C1")