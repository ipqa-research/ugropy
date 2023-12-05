import pytest

from ugropy import instantiate_chem_object


def test_making_it_explode():
    with pytest.raises(ValueError):
        instantiate_chem_object("snake", "mol")

    with pytest.raises(ValueError):
        instantiate_chem_object("acetone", "Argentina")
