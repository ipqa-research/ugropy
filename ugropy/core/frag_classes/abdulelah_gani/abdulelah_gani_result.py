"""Abdulelah-Gani fragmentation model result module :cite:p:`gani`."""

import numpy as np

import pandas as pd

import pint

from rdkit import Chem
from rdkit.Chem import Descriptors

from ugropy.constants import ureg
from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_pst_result import (
    AGaniPSTFragmentationResult,
)


class AGaniFragmentationResult:
    """Abdulelah-Gani fragmentation model result :cite:p:`gani`.

    Class to store the results of the Abdulelah-Gani fragmentation model.

    Parameters
    ----------
    molecule : Chem.rdchem.Mol
        RDKit molecule object.
    primary_fragmentation : AGaniPSTFragmentationResult
        Primary fragmentation model result.
    secondary_fragmentation : AGaniPSTFragmentationResult
        Secondary fragmentation model result.
    tertiary_fragmentation : AGaniPSTFragmentationResult
        Tertiary fragmentation model result.
    properties_contributions : pd.DataFrame
        Contribution parameters of each group for each property of the model.
    properties_biases : pd.DataFrame
        Biases parameters of each property of the model.

    Attributes
    ----------
    molecule : Chem.rdchem.Mol
        RDKit molecule object.
    primary : AGaniPSTFragmentationResult
        Primary fragmentation model result.
    secondary : AGaniPSTFragmentationResult
        Secondary fragmentation model result.
    tertiary : AGaniPSTFragmentationResult
        Tertiary fragmentation model result.
    ml_vector : np.ndarray
        Vector of groups occurrences to evaluate ML model.
    molecular_weight : pint.Quantity
        Molecular weight [g/mol].
    critical_temperature : pint.Quantity
        Critical temperature [K] GC-Simple method.
    critical_pressure : pint.Quantity
        Critical pressure [bar] GC-Simple method.
    critical_volume : pint.Quantity
        Critical volume [cm³/mol] GC-Simple method.
    acentric_factor : pint.Quantity
        Acentric factor [-] GC-Simple method.
    liquid_molar_volume : pint.Quantity
        Liquid molar volume [L/mol] GC-Simple method.
    ig_formation_enthalpy : pint.Quantity
        Ideal gas formation enthalpy [kJ/mol] GC-Simple method.
    ig_formation_gibbs : pint.Quantity
        Ideal gas formation Gibbs [kJ/mol] GC-Simple method.
    """

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
        self.molecular_weight: pint.Quantity = (
            Descriptors.MolWt(self.molecule) * ureg.g / ureg.mol
        )
        self.critical_temperature: pint.Quantity = None
        self.critical_pressure: pint.Quantity = None
        self.critical_volume: pint.Quantity = None
        self.acentric_factor: pint.Quantity = None
        self.liquid_molar_volume: pint.Quantity = None
        self.ig_formation_enthalpy: pint.Quantity = None
        self.ig_formation_gibbs: pint.Quantity = None

        if self.primary.subgroups != {}:
            self.properties_calculation(
                properties_contributions, properties_biases
            )

    def properties_calculation(
        self,
        properties_contributions: pd.DataFrame,
        properties_biases: pd.DataFrame,
    ) -> None:
        """Calculate the properties with the fragmentation results.

        Properties are calculated with the GC-Simple method.

        Parameters
        ----------
        properties_contributions : pd.DataFrame
            Contribution parameters of each group for each property of the
            model.
        properties_biases : pd.DataFrame
            Biases parameters of each property of the model.
        """
        # Critical temperature
        tc_c = properties_contributions["Tc"].values
        tc_b = properties_biases["Tc"].loc["b1"]

        tc_sum = np.dot(tc_c, self.ml_vector.T)[0]

        if tc_sum <= 0:
            self.critical_temperature = np.nan * ureg.K
        else:
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
        w = properties_contributions["w"].values
        w_b1 = properties_biases["w"].loc["b1"]
        w_b2 = properties_biases["w"].loc["b2"]
        w_b3 = properties_biases["w"].loc["b3"]

        w_sum = np.dot(w, self.ml_vector.T)[0]

        self.acentric_factor = (
            np.log((w_sum + w_b3) ** (1 / w_b2)) * w_b1 * ureg.dimensionless
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
