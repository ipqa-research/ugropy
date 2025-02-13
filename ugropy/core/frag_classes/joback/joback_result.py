"""JobackFragmentationResult module."""

from typing import Union

import numpy as np

import pandas as pd

import pint

from rdkit import Chem
from rdkit.Chem import Descriptors

from ugropy.constants import R, ureg
from ugropy.core.frag_classes.base.fragmentation_result import (
    FragmentationResult,
)


class JobackFragmentationResult(FragmentationResult):
    """Joback group contribution properties estimator.

    The class recieves either the Joback and Reid model's :cite:p:`joback1,
    joback2` groups, name or smiles of a molecule and estimates the its
    properties.

    Parameters
    ----------
    molecule : Chem.rdchem.Mol
        RDKit molecule object.
    subgroups : dict
        Dictionary of subgroups.
    subgroups_atoms_indexes : dict
        Dictionary of subgroups atoms indexes.
    properties_contributions : pd.DataFrame
        DataFrame with Joback's properties contributions.
    normal_boiling_point : float, optional
        User provided experimental normal boiling point [K].

    Attributes
    ----------
    subgroups : dict
        Joback functional groups of the molecule.
    experimental_boiling_temperature : pint.Quantity
        User provided experimental normal boiling point [K].
    critical_temperature : pint.Quantity
        Joback estimated critical temperature [K].
    critical_pressure : pint.Quantity
        Joback estimated critical pressure [bar].
    critical_volume : pint.Quantity
        Joback estimated critical volume [cm³/mol].
    normal_boiling_point : pint.Quantity
        Joback estimated normal boiling point [K].
    fusion_temperature : pint.Quantity
        Joback estimated fusion temperature [K].
    ig_formation_formation : pint.Quantity
        Joback estimated enthalpy of formation ideal gas at 298 K [kJ/mol].
    ig_gibbs_formation : pint.Quantity
        Joback estimated Gibbs energy of formation ideal gas at 298 K [K].
    heat_capacity_ideal_gas_params : dict
        Joback estimated Reid's ideal gas heat capacity equation parameters
        [J/mol/K].
    fusion_enthalpy : pint.Quantity
        Joback estimated fusion enthalpy [kJ/mol].
    vaporization_enthalpy : pint.Quantity
        Joback estimated vaporization enthalpy at the normal boiling point
        [kJ/mol].
    sum_na : float
        Joback n_A contribution to liquid viscosity [Pa s].
    sum_nb : float
        Joback n_B contribution to liquid viscosity [Pa s].
    molecular_weight : pint.Quantity
        Molecular weight from Joback's subgroups [g/mol].
    acentric_factor : pint.Quantity
        Acentric factor from Lee and Kesler's equation :cite:p:`joback1`.
    vapor_pressure_params : dict
        Vapor pressure G and k parameters for the Riedel-Plank-Miller
        equation [bar] :cite:p:`joback1`.
    """

    def __init__(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
        properties_contributions: pd.DataFrame,
        normal_boiling_point: float = None,
    ) -> None:
        # Initialize the FragmentationResult attributes
        super().__init__(molecule, subgroups, subgroups_atoms_indexes)

        # experimental boiling temperature
        self.experimental_boiling_temperature = normal_boiling_point

        # Original Joback properties
        self.critical_temperature: pint.Quantity = None
        self.critical_pressure: pint.Quantity = None
        self.critical_volume: pint.Quantity = None
        self.normal_boiling_point: pint.Quantity = None
        self.fusion_temperature: pint.Quantity = None
        self.ig_enthalpy_formation: pint.Quantity = None
        self.ig_gibbs_formation: pint.Quantity = None
        self.heat_capacity_ideal_gas_params = np.array([])
        self.fusion_enthalpy: pint.Quantity = None
        self.vaporization_enthalpy: pint.Quantity = None
        self.sum_na = None
        self.sum_nb = None
        self.molecular_weight: pint.Quantity = None

        # Extra properties
        self.acentric_factor: pint.Quantity = None
        self.vapor_pressure_params = {}

        # Fill the properties' values
        if self.subgroups != {}:
            self._calculate_properties(properties_contributions)

    def heat_capacity_ideal_gas(
        self, temperature: Union[float, np.ndarray]
    ) -> pint.Quantity:
        """Calculate the ideal gas heat capacity [J/mol/K].

        Uses the Joback estimated Reid's ideal gas heat capacity equation
        parameters [J/mol/K].

        Parameters
        ----------
        temperature : Union[float, np.ndarray]
            Temperature [K]

        Returns
        -------
        Union[float, np.ndarray]
            Ideal gas heat capacity [J/mol/K].
        """
        a, b, c, d = self.heat_capacity_ideal_gas_params

        t = temperature

        return (a + b * t + c * t**2 + d * t**3) * ureg.J / ureg.mol / ureg.K

    def heat_capacity_liquid(
        self, temperature: Union[float, np.ndarray]
    ) -> pint.Quantity:
        """Calculate the liquid heat capacity [J/mol/K].

        Uses the Rowlinson-Bondi :cite:p:`joback1` equation with the Joback
        estimated properties.

        Parameters
        ----------
        temperature : Union[float, np.ndarray]
            Temperature [K]

        Returns
        -------
        Union[float, np.ndarray]
            Ideal gas heat capacity [J/mol/K].
        """
        tr = temperature / self.critical_temperature.magnitude
        w = self.acentric_factor.magnitude

        c_p0 = self.heat_capacity_ideal_gas(temperature).magnitude

        c_pl = c_p0 + R * (
            2.56
            + 0.436 * (1 - tr) ** (-1)
            + w
            * (
                2.91
                + 4.28 * (1 - tr) ** (-1 / 3) * tr ** (-1)
                + 0.296 * (1 - tr) ** (-1)
            )
        )

        return c_pl * ureg.J / ureg.mol / ureg.K

    def viscosity_liquid(
        self, temperature: Union[float, np.ndarray]
    ) -> pint.Quantity:
        """Calculate the Joback estimated liquid viscosity [Pa s].

        Parameters
        ----------
        temperature : Union[float, np.ndarray]
            Temperature [K]

        Returns
        -------
        Union[float, np.ndarray]
            Liquid viscosity [Pa s].
        """
        t = temperature

        n_l = self.molecular_weight.magnitude * np.exp(
            (self.sum_na - 597.82) / t + self.sum_nb - 11.202
        )
        return n_l * ureg.Pa * ureg.s

    def vapor_pressure(
        self, temperature: Union[float, np.ndarray]
    ) -> pint.Quantity:
        """Calculate the vapor pressure [bar].

        Uses the Riedel-Plank-Miller :cite:p:`joback1` equation with the Joback
        estimated properties.

        Parameters
        ----------
        temperature : Union[float, np.ndarray]
            Temperature [K]

        Returns
        -------
        Union[float, np.ndarray]
            Vapor pressure [bar]
        """
        tr = temperature / self.critical_temperature.magnitude

        g = self.vapor_pressure_params["G"]
        k = self.vapor_pressure_params["k"]

        vp_r = np.exp(-g / tr * (1 - tr**2 + k * (3 + tr) * (1 - tr) ** 3))

        vp = vp_r * self.critical_pressure.magnitude

        return vp * ureg.bar

    def _calculate_properties(self, contribs: pd.DataFrame) -> None:
        """Obtain the molecule properties from Joback's groups."""
        groups = list(self.subgroups.keys())
        ocurr = list(self.subgroups.values())

        df = contribs.loc[groups]

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
        mw_c = Descriptors.MolWt(self.molecule)

        # Molecular weight
        self.molecular_weight = mw_c * ureg.g / ureg.mol

        # Joback normal boiling point (Tb)
        self.normal_boiling_point = (198.2 + np.dot(ocurr, tb_c)) * ureg.K

        # Fusion temperature (Tf)
        self.fusion_temperature = (122.5 + np.dot(ocurr, tf_c)) * ureg.K

        # Used normal boiling point for calculations
        if self.experimental_boiling_temperature is not None:
            tb = self.experimental_boiling_temperature
        else:
            tb = self.normal_boiling_point.magnitude

        # Critical temperature (Tc) normal boiling temperature for calculations
        self.critical_temperature = (
            tb
            * (
                0.584
                + 0.965 * np.dot(ocurr, tc_c)
                - (np.dot(ocurr, tc_c)) ** 2
            )
            ** (-1)
            * ureg.K
        )

        # Critical pressure (Pc)
        self.critical_pressure = (
            0.113 + 0.0032 * np.dot(ocurr, numa_c) - np.dot(ocurr, pc_c)
        ) ** (-2) * ureg.bar

        # Critical volume (Vc)
        self.critical_volume = (
            (17.5 + np.dot(ocurr, vc_c)) * ureg.cm**3 / ureg.mol
        )

        # Standard enthalpy of formation (298 K)
        self.ig_enthalpy_formation = (
            (68.29 + np.dot(ocurr, hform_c)) * ureg.kJ / ureg.mol
        )

        # Standard Gibbs energy of formation (298 K)
        self.ig_gibbs_formation = (
            (53.88 + np.dot(ocurr, gform_c)) * ureg.kJ / ureg.mol
        )

        # Enthalpy of vaporization
        self.vaporization_enthalpy = (
            (15.30 + np.dot(ocurr, hvap_c)) * ureg.kJ / ureg.mol
        )

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

            self.heat_capacity_ideal_gas_params = np.array([a, b, c, d])

        # Enthalpy of fusion
        if all(df["Hfusion"].notnull()):
            hfusion_c = df["Hfusion"].to_numpy()
            self.fusion_enthalpy = (
                (-0.88 + np.dot(ocurr, hfusion_c)) * ureg.kJ / ureg.mol
            )

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
        t_br = tb / self.critical_temperature.magnitude

        # Lee and Kesler's equation (acentric factor)
        pc = self.critical_pressure.magnitude

        self.acentric_factor = (
            (
                -np.log(pc)
                - 5.92714
                + 6.09648 / t_br
                + 1.28862 * np.log(t_br)
                - 0.169347 * t_br**6
            )
            / (
                15.2518
                - 15.6875 / t_br
                - 13.4721 * np.log(t_br)
                + 0.43577 * t_br**6
            )
            * ureg.dimensionless
        )

        # Riedel-Plank-Miller equation (vapor pressure [bar])
        h = t_br * np.log(pc / 1.01325) / (1 - t_br)

        g = 0.4835 + 0.4605 * h

        k = (h / g - (1 + t_br)) / ((3 + t_br) * (1 - t_br) ** 2)

        self.vapor_pressure_params = {"G": g, "k": k}
