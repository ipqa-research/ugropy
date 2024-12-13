"""Groups module."""

from typing import List, Union

from rdkit.Chem import Descriptors

from .core.frag_classes.gibbs_model.gibbs_result import (
    GibbsFragmentationResult,
)
from .core.frag_classes.joback.joback_result import JobackFragmentationResult
from .core.get_rdkit_object import instantiate_mol_object
from .core.ilp_solvers.default_solver import DefaultSolver
from .core.ilp_solvers.ilp_solver import ILPSolver
from .models.jobackmod import joback
from .models.psrkmod import psrk
from .models.unifacmod import unifac


class Groups:
    """Group class.

    Stores the solved FragmentationModels subgroups of a molecule.

    Parameters
    ----------
    identifier : str or rdkit.Chem.rdchem.Mol
        Identifier of a molecule (name, SMILES or Chem.rdchem.Mol). Example:
        hexane or CCCCCC.
    identifier_type : str, optional
        Use 'name' to search a molecule by name, 'smiles' to provide the
        molecule SMILES representation or 'mol' to provide a
        rdkit.Chem.rdchem.Mol object, by default "name".
    solver : ILPSolver, optional
        ILP solver to use, by default DefaultSolver.
    search_multiple_solutions : bool, optional
        If True, the solver will search for multiple solutions. If set to true,
        the model's results will be lists of FragmentationResults objects, by
        default False.
    normal_boiling_temperature : float, optional
        If provided, will be used to estimate critical temperature, acentric
        factor, and vapor pressure instead of the estimated normal boiling
        point in the Joback group contribution model, by default None.

    Attributes
    ----------
    identifier : str
        Identifier of a molecule. Example: hexane or CCCCCC.
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
    joback : Union[JobackFragmentationResult, List[JobackFragmentationResult]]
        JobackFragmentationResult object that contains the Joback subgroups and
        the estimated properties of the molecule.
    """

    def __init__(
        self,
        identifier: str,
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
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
        )

        # PSRK
        self.psrk: Union[
            GibbsFragmentationResult, List[GibbsFragmentationResult]
        ] = psrk.get_groups(
            self.identifier,
            self.identifier_type,
            solver=solver,
            search_multiple_solutions=search_multiple_solutions,
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
        )
