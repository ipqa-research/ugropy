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
        normal_boiling_temperature: float = None,
    ) -> None:
        # Skip if instantiation from_groups is made.
        if identifier != "__skip__":
            self.groups = get_joback_groups(identifier, identifier_type)
            self.exp_nbt = normal_boiling_temperature
        else:
            self.groups = {}
            self.exp_nbt = None

        # Original Joback properties
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

        # Extra properties
        self.acentric_factor = None
        self.vapor_pressure_params = {}

        # Fill the properties' values
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

    def vapor_pressure(self, temperature):
        tr = temperature / self.critical_temperature

        g = self.vapor_pressure_params["G"]
        k = self.vapor_pressure_params["k"]

        vp_r = np.exp(-g / tr * (1 - tr**2 + k * (3 + tr) * (1 - tr)**3))

        vp = vp_r * self.critical_pressure

        return vp

    @classmethod
    def from_groups(
        cls, joback_groups: dict, normal_boiling_temperature: float = None
    ) -> "Joback":
        mol = cls("__skip__")
        mol.groups = joback_groups
        mol.exp_nbt = normal_boiling_temperature
        mol._calculate_properties()
        return mol

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
        # normal boiling temperature for calculations
        if self.exp_nbt is not None:
            tb = self.exp_nbt
        else:
            tb = self.normal_boiling_temperature

        self.critical_temperature = tb * (
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

        # =====================================================================
        # Extra properties
        # =====================================================================
        # Reduced normal boiling point temperature
        t_br = tb / self.critical_temperature

        # Lee and Kesler's equation (acentric factor)
        pc = self.critical_pressure
        self.acentric_factor = (
            -np.log(pc)
            - 5.92714
            + 6.09648 / t_br
            + 1.28862 * np.log(t_br)
            - 0.169347 * t_br**6
        ) / (
            15.2518
            - 15.6875 / t_br
            - 13.4721 * np.log(t_br)
            + 0.43577 * t_br**6
        )

        # Riedel-Plank-Miller equation (vapor pressure [bar])
        h = t_br * np.log(self.critical_pressure / 1.01325) / (1 - t_br)

        g = 0.4835 + 0.4605 * h

        k = (h / g - (1 + t_br)) / ((3 + t_br) * (1 - t_br)**2)

        self.vapor_pressure_params = {"G": g, "k": k}
