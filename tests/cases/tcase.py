from abc import ABC, abstractmethod

import pytest

from ugropy import DefaultSolver, ILPSolver

from .acids import acids_cases
from .alcohols import alcohols_cases
from .aldehydes import aldehydes_cases
from .aromatics import aromatics_cases
from .case import Case
from .complex import complex_cases
from .epoxides import epoxides_cases
from .esters import esters_cases
from .ethers import ethers_cases
from .halogens import halogens_cases
from .hydrocarbons import hydrocarbons_cases
from .ketones import ketones_cases
from .nitrogen import nitrogen_cases
from .particulars import particulars_cases
from .silicon import silicon_cases
from .sulfur import sulfur_cases
from .unsaturated_hydrocarbons import unsaturated_hydrocarbons_cases


# Having multiple solvers in the future we can check them all here.
solvers = [DefaultSolver]


# Class to define the test class for each model.
class TCase(ABC):
    """Class to define the test class for each model.

    This class is an abstract class that defines the structure of the test
    class for each model. It defines the asserts method that should be
    implemented in the subclasses, and the test methods for each type of
    molecule. The test methods are parametrized with the solvers and the cases
    of each type of molecule.
    """

    @abstractmethod
    def asserts(self, case: Case, solver: ILPSolver) -> None: ...

    @pytest.mark.acids
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", acids_cases, ids=[c.identifier for c in acids_cases]
    )
    def test_acids(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.alcohols
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", alcohols_cases, ids=[c.identifier for c in alcohols_cases]
    )
    def test_alcohols(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.aldehydes
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", aldehydes_cases, ids=[c.identifier for c in aldehydes_cases]
    )
    def test_aldehydes(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.aromatics
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", aromatics_cases, ids=[c.identifier for c in aromatics_cases]
    )
    def test_aromathics(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.complex
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", complex_cases, ids=[c.identifier for c in complex_cases]
    )
    def test_complex(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.epoxides
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", epoxides_cases, ids=[c.identifier for c in epoxides_cases]
    )
    def test_epoxides(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.esters
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", esters_cases, ids=[c.identifier for c in esters_cases]
    )
    def test_esters(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.ethers
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", ethers_cases, ids=[c.identifier for c in ethers_cases]
    )
    def test_ethers(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.halogens
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", halogens_cases, ids=[c.identifier for c in halogens_cases]
    )
    def test_halogens(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.hydrocarbons
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case",
        hydrocarbons_cases,
        ids=[c.identifier for c in hydrocarbons_cases],
    )
    def test_hydrocarbons(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.unsaturated_hydrocarbons
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case",
        unsaturated_hydrocarbons_cases,
        ids=[c.identifier for c in unsaturated_hydrocarbons_cases],
    )
    def test_insaturated_hydrocarbons(
        self, case: Case, solver: ILPSolver
    ) -> None:
        self.asserts(case, solver)

    @pytest.mark.ketones
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", ketones_cases, ids=[c.identifier for c in ketones_cases]
    )
    def test_ketones(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.nitrogen
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", nitrogen_cases, ids=[c.identifier for c in nitrogen_cases]
    )
    def test_nitrogen(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.particulars
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case",
        particulars_cases,
        ids=[c.identifier for c in particulars_cases],
    )
    def test_particulars(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.silicon
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", silicon_cases, ids=[c.identifier for c in silicon_cases]
    )
    def test_silicon(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)

    @pytest.mark.sulfur
    @pytest.mark.parametrize("solver", solvers)
    @pytest.mark.parametrize(
        "case", sulfur_cases, ids=[c.identifier for c in sulfur_cases]
    )
    def test_sulfur(self, case: Case, solver: ILPSolver) -> None:
        self.asserts(case, solver)
