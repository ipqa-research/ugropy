import numpy as np

from ugropy.constants import joback_properties_contibutions, joback_subgroups
from ugropy.core.model_getters import (
    get_joback_groups,
)


class Joback:
    def __init__(
        self,
        identifier: str,
        identifier_type: str = "name",
    ) -> None:
        self.groups = get_joback_groups(identifier, identifier_type)

        self.critical_temperature = None
        self.critical_pressure = None
        self.critical_volume = None
        self.normal_boiling_temperature = None
        self.fusion_temperature = None
        self.h_formation = None
        self.g_formation = None
        self.heat_capacity_params = None
        self.h_fusion = None
        self.h_vaporization = None
        self.sum_na = None
        self.sum_nb = None
        self.molecular_weight = None

        if self.groups != {}:
            self._calculate_properties()

    def heat_capacity(self, temperature):
        a, b, c, d = self.heat_capacity_params

        t = temperature

        return a + b * t + c * t**2 + d * t**3

    def liquid_viscosity(self, temperature):
        t = temperature

        n_l = self.molecular_weight * np.exp(
            (self.sum_na - 597.82) / t + self.sum_nb - 11.202
        )
        return n_l

    def _calculate_properties(self):
        groups = list(self.groups.keys())
        ocurr = list(self.groups.values())

        df = joback_properties_contibutions.loc[groups]

        # =====================================================================
        # Calculate complete contribution properties (no contribution missing)
        # =====================================================================
        tc_c = df["Tc"].to_numpy()
        pc_c = df["Pc"].to_numpy()
        vc_c = df["Vc"].to_numpy()
        tb_c = df["Tb"].to_numpy()
        tf_c = df["Tf"].to_numpy()
        hform_c = df["Hform"].to_numpy()
        gform_c = df["Gform"].to_numpy()
        hvap_c = df["Hvap"].to_numpy()
        numa_c = df["num_of_atoms"].to_numpy()
        mw_c = joback_subgroups.loc[groups, "molecular_weight"].to_numpy()

        # Molecular weight
        self.molecular_weight = np.dot(ocurr, mw_c)

        # Joback normal boiling point (Tb)
        self.normal_boiling_temperature = 198.2 + np.dot(ocurr, tb_c)

        # Fusion temperature (Tf)
        self.fusion_temperature = 122.5 + np.dot(ocurr, tf_c)

        # Critical temperature (Tc)
        self.critical_temperature = self.normal_boiling_temperature * (
            0.584 + 0.965 * np.dot(ocurr, tc_c) - (np.dot(ocurr, tc_c)) ** 2
        ) ** (-1)

        # Critical pressure (Pc)
        self.critical_pressure = (
            0.113 + 0.0032 * np.dot(ocurr, numa_c) - np.dot(ocurr, pc_c)
        ) ** (-2)

        # Critical volume (Vc)
        self.critical_volume = 17.5 + np.dot(ocurr, vc_c)

        # Standard enthalpy of formation (298 K)
        self.h_formation = 68.29 + np.dot(ocurr, hform_c)

        # Standard Gibbs energy of formation (298 K)
        self.g_formation = 53.88 + np.dot(ocurr, gform_c)

        # Enthalpy of vaporization
        self.h_vaporization = 15.30 + np.dot(ocurr, hvap_c)

        # =====================================================================
        # Incomplete contribution properties (some contribution missing)
        # =====================================================================
        # Heat capacity
        if "-N= (ring)" not in groups:
            a_c = df["a"].to_numpy()
            b_c = df["b"].to_numpy()
            c_c = df["c"].to_numpy()
            d_c = df["d"].to_numpy()

            a = np.dot(ocurr, a_c) - 37.93
            b = np.dot(ocurr, b_c) + 0.21
            c = np.dot(ocurr, c_c) - 3.91e-4
            d = np.dot(ocurr, d_c) + 2.06e-7

            self.heat_capacity_params = np.array([a, b, c, d])

        # Enthalpy of fusion
        if all(df["Hfusion"].notnull()):
            hfusion_c = df["Hfusion"].to_numpy()
            self.h_fusion = -0.88 + np.dot(ocurr, hfusion_c)

        # Liquid viscosity
        if all(df["na"].notnull()):
            na_c = df["na"].to_numpy()
            nb_c = df["nb"].to_numpy()

            self.sum_na = np.dot(ocurr, na_c)
            self.sum_nb = np.dot(ocurr, nb_c)
