import numpy as np

import pandas as pd

import pint

from rdkit import Chem

from ugropy.constants import ureg
from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_pst_result import (
    AGaniPSTFragmentationResult,
)


class AGaniFragmentationResult:
    def __init__(
        self,
        molecule: Chem.rdchem.Mol,
        primary_fragmentation: AGaniPSTFragmentationResult,
        secondary_fragmentation: AGaniPSTFragmentationResult,
        tertiary_fragmentation: AGaniPSTFragmentationResult,
        properties_contributions: pd.DataFrame,
        properties_biases: pd.DataFrame,
    ) -> None:

        # =====================================================================
        # Fragmentation models
        # =====================================================================
        self.molecule = molecule
        self.primary = primary_fragmentation
        self.secondary = secondary_fragmentation
        self.tertiary = tertiary_fragmentation

        # =====================================================================
        # ML Vector. Vector of parameters occurrences to evaluate ML model
        # =====================================================================
        self.ml_vector = np.zeros(424, dtype=np.int64)

        for group, ocurrences in self.primary.subgroups_numbers.items():
            self.ml_vector[group - 1] = ocurrences

        for group, ocurrences in self.secondary.subgroups_numbers.items():
            self.ml_vector[group - 1] = ocurrences

        for group, ocurrences in self.tertiary.subgroups_numbers.items():
            self.ml_vector[group - 1] = ocurrences

        self.ml_vector = self.ml_vector.reshape(1, -1)

        # =====================================================================
        # Properties
        # =====================================================================
        self.critical_temperature: pint.Quantity = None
        self.critical_pressure: pint.Quantity = None
        self.critical_volume: pint.Quantity = None
        self.acentric_factor: pint.Quantity = None
        self.liquid_molar_volume: pint.Quantity = None
        self.ig_formation_enthalpy: pint.Quantity = None
        self.ig_formation_gibbs: pint.Quantity = None
        self.normal_melting_point: pint.Quantity = None
        self.normal_boiling_point: pint.Quantity = None

        if self.primary.subgroups != {}:
            self.properties_calculation(
                properties_contributions, properties_biases
            )

    def properties_calculation(
        self, properties_contributions, properties_biases
    ) -> None:
        # Critical temperature
        tc_c = properties_contributions["Tc"].values
        tc_b = properties_biases["Tc"].loc["b1"]

        tc_sum = np.dot(tc_c, self.ml_vector.T)[0]

        self.critical_temperature = tc_b * np.log(tc_sum) * ureg.K

        # Critical pressure
        pc_c = properties_contributions["Pc"].values
        pc_b1 = properties_biases["Pc"].loc["b1"]
        pc_b2 = properties_biases["Pc"].loc["b2"]

        pc_sum = np.dot(pc_c, self.ml_vector.T)[0]

        self.critical_pressure = (
            (pc_sum + pc_b2) ** (-1 / 0.5) + pc_b1
        ) * ureg.bar

        # Critical volume
        vc_c = properties_contributions["Vc"].values
        vc_b1 = properties_biases["Vc"].loc["b1"]

        vc_sum = np.dot(vc_c, self.ml_vector.T)[0]

        self.critical_volume = (vc_sum + vc_b1) * ureg.cm**3 / ureg.mol

        # Acentric factor
        w = properties_contributions["ω"].values
        w_b1 = properties_biases["ω"].loc["b1"]
        w_b2 = properties_biases["ω"].loc["b2"]
        w_b3 = properties_biases["ω"].loc["b3"]

        w_sum = np.dot(w, self.ml_vector.T)[0]

        self.acentric_factor = (
            np.log((w_sum + w_b3)) ** (1 / w_b2) * w_b1 * ureg.dimensionless
        )

        # Liquid molar volume
        lmv_c = properties_contributions["Lmv"].values
        lmv_b1 = properties_biases["Lmv"].loc["b1"]

        lmv_sum = np.dot(lmv_c, self.ml_vector.T)[0]

        self.liquid_molar_volume = (lmv_sum + lmv_b1) * ureg.L / ureg.mol

        # Ideal gas formation enthalpy
        h_f_c = properties_contributions["Hf"].values
        h_f_c_b1 = properties_biases["Hf"].loc["b1"]

        h_f_sum = np.dot(h_f_c, self.ml_vector.T)[0]

        self.ig_formation_enthalpy = (h_f_sum + h_f_c_b1) * ureg.kJ / ureg.mol

        # Ideal gas formation Gibbs
        g_f_c = properties_contributions["Gf"].values
        g_f_c_b1 = properties_biases["Gf"].loc["b1"]

        g_f_sum = np.dot(g_f_c, self.ml_vector.T)[0]

        self.ig_formation_gibbs = (g_f_sum + g_f_c_b1) * ureg.kJ / ureg.mol

        # Normal melting point
        tm_c = properties_contributions["Tm"].values
        tm_b1 = properties_biases["Tm"].loc["b1"]

        tm_sum = np.dot(tm_c, self.ml_vector.T)[0]

        self.normal_melting_point = np.log(tm_sum) * tm_b1 * ureg.K

        # Normal boiling point
        tb_c = properties_contributions["Tb"].values
        tb_b1 = properties_biases["Tb"].loc["b1"]

        tb_sum = np.dot(tb_c, self.ml_vector.T)[0]

        self.normal_boiling_point = np.log(tb_sum) * tb_b1 * ureg.K
