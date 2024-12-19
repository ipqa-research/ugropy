import numpy as np

from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_pst_result import (
    AGaniPSTFragmentationResult,
)


class AGaniFragmentationResult:
    def __init__(
        self,
        primary_fragmentation: AGaniPSTFragmentationResult,
        secondary_fragmentation: AGaniPSTFragmentationResult,
        tertiary_fragmentation: AGaniPSTFragmentationResult,
    ) -> None:

        self.primary = primary_fragmentation
        self.secondary = secondary_fragmentation
        self.tertiary = tertiary_fragmentation

        self.groups_vector = np.zeros(424, dtype=np.int64)

        for group, ocurrences in self.primary.subgroups_numbers.items():
            self.groups_vector[group - 1] = ocurrences

        for group, ocurrences in self.secondary.subgroups_numbers.items():
            self.groups_vector[group - 1] = ocurrences

        for group, ocurrences in self.tertiary.subgroups_numbers.items():
            self.groups_vector[group - 1] = ocurrences
            
        self.groups_vector = self.groups_vector.reshape(1, -1)
