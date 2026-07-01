"""Groups class module."""

from typing import List, Union

from rdkit.Chem import Descriptors

from .core.frag_classes.abdulelah_gani.abdulelah_gani_result import (
    AGaniFragmentationResult,
)
from .core.frag_classes.gibbs_model.gibbs_result import (
    GibbsFragmentationResult,
)
from .core.frag_classes.joback.joback_result import JobackFragmentationResult
from .core.get_rdkit_object import instantiate_mol_object
from .core.ilp_solvers.default_solver import DefaultSolver
from .core.ilp_solvers.ilp_solver import ILPSolver
from .models.abdulelah_gani_mod import abdulelah_gani
from .models.dortmundmod import dortmund
from .models.jobackmod import joback
from .models.psrkmod import psrk
from .models.unifacmod import unifac


class Groups:
    """Groups class.

    Stores the solved FragmentationModels results of a molecule. This class
    was implemented on an early version of the `ugropy` library. Is not really
    meant to be used, instead is recommended to use directly the corresponding
    FragmentationModel independently since it provides more flexibility and
    control over the results. The class is kept since in most cases is very
    comfortable to have all the results in a single object.

    Parameters
    ----------
    identifier : Union[str, Chem.rdchem.Mol]
        Identifier of the molecule. You can use either the name of the
        molecule, the SMILEs of the molecule or a rdkit Mol object.
    identifier_type : str, optional
        Identifier type of the molecule. Use "name" if you are providing
        the molecules' name, "smiles" if you are providing the SMILES
        or "mol" if you are providing a rdkir mol object, by default "name"
    solver : ILPSolver, optional
        ILP solver class, by default DefaultSolver
    search_multiple_solutions : bool, optional
        Weather search for multiple solutions or not, by default False
        If False the return will be a FragmentationResult object, if True
        the return will be a list of FragmentationResult objects.
    search_nonoptimal : bool, optional
        If True, the solver will search for non-optimal solutions along with
        the optimal ones. This is useful when the user wants to find all
        possible combinations of fragments that cover the universe. By default
        False. If `search_multiple_solutions` is False, this parameter will be
        ignored.
    normal_boiling_temperature : float, optional
        If provided, will be used to estimate critical temperature, acentric
        factor, and vapor pressure instead of the estimated normal boiling
        point in the Joback group contribution model, by default None.

    Attributes
    ----------
    identifier : str
        Identifier of a molecule. Example: 'hexane' or 'CCCCCC'.
    identifier_type : str, optional
        Use 'name' to search a molecule by name or 'smiles' to provide the
        molecule SMILES representation, by default "name".
    mol_object : rdkit.Chem.rdchem.Mol
        RDKit Mol object.
    molecular_weight : float
        Molecule's molecular weight from rdkit.Chem.Descriptors.MolWt [g/mol].
    unifac : Union[GibbsFragmentationResult, List[GibbsFragmentationResult]]
        Classic LV-UNIFAC subgroups.
    psrk : Union[GibbsFragmentationResult, List[GibbsFragmentationResult]]
        Predictive Soave-Redlich-Kwong subgroups.
    dortmund : Union[GibbsFragmentationResult, List[GibbsFragmentationResult]]
        Dortmund UNIFAC subgroups.
    joback : Union[JobackFragmentationResult, List[JobackFragmentationResult]]
        JobackFragmentationResult object that contains the Joback subgroups and
        the estimated properties of the molecule.
    agani : Union[AGaniFragmentationResult, List[AGaniFragmentationResult]]
        AGaniFragmentationResult object that contains the Abdulelah-Gani
        subgroups and the estimated properties of the molecule.
    """

    def __init__(
        self,
        identifier: str,
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
        search_nonoptimal: bool = False,
        normal_boiling_temperature: float = None,
    ) -> None:
        self.identifier_type = identifier_type.lower()
        self.identifier = identifier
        self.mol_object = instantiate_mol_object(identifier, identifier_type)
        self.molecular_weight = Descriptors.MolWt(self.mol_object)

        # UNIFAC
        self.unifac: Union[
            GibbsFragmentationResult, List[GibbsFragmentationResult]
        ] = unifac.get_groups(
            self.identifier,
            self.identifier_type,
            solver=solver,
            search_multiple_solutions=search_multiple_solutions,
            search_nonoptimal=search_nonoptimal,
        )

        # PSRK
        self.psrk: Union[
            GibbsFragmentationResult, List[GibbsFragmentationResult]
        ] = psrk.get_groups(
            self.identifier,
            self.identifier_type,
            solver=solver,
            search_multiple_solutions=search_multiple_solutions,
            search_nonoptimal=search_nonoptimal,
        )

        # Dortmund
        self.dortmund: Union[
            GibbsFragmentationResult, List[GibbsFragmentationResult]
        ] = dortmund.get_groups(
            self.identifier,
            self.identifier_type,
            solver=solver,
            search_multiple_solutions=search_multiple_solutions,
            search_nonoptimal=search_nonoptimal,
        )

        # Joback
        self.joback: Union[
            JobackFragmentationResult, List[JobackFragmentationResult]
        ] = joback.get_groups(
            self.identifier,
            self.identifier_type,
            solver=solver,
            search_multiple_solutions=search_multiple_solutions,
            normal_boiling_point=normal_boiling_temperature,
            search_nonoptimal=search_nonoptimal,
        )

        # Abdulelah-Gani
        self.agani: Union[
            AGaniFragmentationResult, List[AGaniFragmentationResult]
        ] = abdulelah_gani.get_groups(
            self.identifier,
            self.identifier_type,
            solver=solver,
            search_multiple_solutions=search_multiple_solutions,
            search_nonoptimal=search_nonoptimal,
        )
