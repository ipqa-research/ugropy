from .get_groups import get_groups

from .constants import unifac_subgroups, unifac_matrix

class Substance:
    def __init__(self, name) -> None:
        self.name = name
        self.unifac_groups = get_groups(self.name, unifac_subgroups, unifac_matrix)