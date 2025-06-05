"""Joback fragmentation module."""

from typing import List, Union

import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_pst_result import (
    AGaniPSTFragmentationResult,
)
from ugropy.core.frag_classes.base.fragmentation_model import (
    FragmentationModel,
)
from ugropy.core.ilp_solvers.default_solver import DefaultSolver
from ugropy.core.ilp_solvers.ilp_solver import ILPSolver


class AbdulelahGaniPSTModel(FragmentationModel):
    """Abdulelah-Gani model dedicated to properties estimation.

    Class to construct the primary, secondary and tertiary structures detector
    for the Abdulelah-Gani properties estimation model :cite:p:`gani`.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule).
    subgroups_info : pd.DataFrame
        Group's subgroups numbers.

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule).
    detection_mols : dict
        Dictionary cotaining all the rdkit Mol object from the detection_smarts
        subgroups column
    info : pd.DataFrame
        Group's subgroups numbers.
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        subgroups_info: pd.DataFrame,
        allow_overlapping: bool = False,
        allow_free_atoms: bool = False,
    ) -> None:

        super().__init__(
            subgroups=subgroups,
            allow_overlapping=allow_overlapping,
            allow_free_atoms=allow_free_atoms,
            fragmentation_result=AGaniPSTFragmentationResult,
        )

        self.subgroups_info = subgroups_info

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
        search_nonoptimal: bool = False,
        solver_arguments: dict = {},
    ) -> Union[AGaniPSTFragmentationResult, List[AGaniPSTFragmentationResult]]:
        """Get the groups of a molecule.

        Parameters
        ----------
        identifier : Union[str, Chem.rdchem.Mol]
            Identifier of the molecule. You can use either the name of the
            molecule, the SMILEs of the molecule or a rdkit Mol object.
        identifier_type : str, optional
            Identifier type of the molecule. Use "name" if you are providing
            the molecules' name, "smiles" if you are providing the SMILES or
            "mol" if you are providing a rdkir mol object, by default "name"
        solver : ILPSolver, optional
            ILP solver class, by default DefaultSolver
        search_multiple_solutions : bool, optional
            Weather search for multiple solutions or not, by default False If
            False the return will be a FragmentationResult object, if True the
            return will be a list of FragmentationResult objects.
        search_nonoptimal : bool, optional
            If True, the solver will search for non-optimal solutions along
            with the optimal ones. This is useful when the user wants to find
            all possible combinations of fragments that cover the universe. By
            default False. If `search_multiple_solutions` is False, this
            parameter will be ignored.
        solver_arguments : dict, optional
            Dictionary with the arguments to be passed to the solver. For the
            DefaultSolver of ugropy you can change de PulP solver passing a
            dictionary like {"solver": "PULP_CBC_CMD"} and change the PulP
            solver. If empty it will use the default solver arguments, by
            default {}.

        Returns
        -------
        Union[AGaniPFragmentationResult, List[AGaniPFragmentationResult]]
            Fragmentation result. If search_multiple_solutions is False the
            return will be a FragmentationResult object, if True the return
            will be a list of FragmentationResult objects.
        """
        sol = super().get_groups(
            identifier,
            identifier_type,
            solver,
            search_multiple_solutions,
            search_nonoptimal=search_nonoptimal,
            solver_arguments=solver_arguments,
            subgroups_info=self.subgroups_info,
        )

        return sol

    def mol_preprocess(self, mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
        """Preprocess the molecule to be ready for the fragmentation.

        This method preprocess the molecule to be ready for the fragmentation
        process. The preprocessing steps are:

        1. Kekulize the molecule.
        2. Identify the aromatic rings.
        3. Check the aromaticity of the rings.
        4. Make the rings aromatic or non-aromatic based on the setp 3 check.

        The criteria to check the aromaticity of the rings are based on the
        criteria proposed by Abdulelah-Gani in the original paper datase
        :cite:p:`gani`.

        Parameters
        ----------
        mol : Chem.rdchem.Mol
            Molecule to preprocess

        Returns
        -------
        Chem.rdchem.Mol
            Preprocessed molecule
        """
        # Clone of the molecule to kekulize
        kekulized_mol = Chem.Mol(mol)
        Chem.Kekulize(kekulized_mol, clearAromaticFlags=True)

        # Step 1: Identify all ring on the molecule
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()

        # Step 2: Get only the aromatic rings
        rdkit_aromatic_rings = [
            ring
            for ring in atom_rings
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]

        # Lists to store the AGani aromatic rings and not aromatic rings
        aromatic_rings = []
        non_aromatic_rings = []

        # Step 3: Check which rings are aromatic
        for ring in rdkit_aromatic_rings:
            if self.check_ring_aromaticity(mol, ring):
                aromatic_rings.append(ring)
            else:
                non_aromatic_rings.append(ring)

        # Step 4: Preprocess the non aromatic rings
        for ring in non_aromatic_rings:
            # Make atoms of the ring non-aromatic
            for idx in ring:
                mol.GetAtomWithIdx(idx).SetIsAromatic(False)

            # Make bonds of the ring non-aromatic
            for i in range(len(ring)):
                atom1, atom2 = ring[i], ring[(i + 1) % len(ring)]
                bond = mol.GetBondBetweenAtoms(atom1, atom2)

                # Set bond type to the kekulized bond type
                kekulized_bond = kekulized_mol.GetBondBetweenAtoms(
                    atom1, atom2
                )
                bond.SetBondType(kekulized_bond.GetBondType())

        # Step 5: Revisit aromatic rings to make sure they are aromatic
        for ring in aromatic_rings:
            # Restore atoms aromaticity
            for idx in ring:
                mol.GetAtomWithIdx(idx).SetIsAromatic(True)

            # Restore bonds aromaticity
            for i in range(len(ring)):
                atom1, atom2 = ring[i], ring[(i + 1) % len(ring)]
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                bond.SetBondType(Chem.rdchem.BondType.AROMATIC)

        return mol

    def check_ring_aromaticity(
        self, mol: Chem.rdchem.Mol, ring: List[int]
    ) -> bool:
        """Verify ring aromaticity along with the authors criteria.

        Parameters
        ----------
        mol : Chem.Mol
            Molecule
        ring : list[int]
            List of atom indexes that form the ring

        Returns
        -------
        bool
            True if the ring is aromatic, False otherwise
        """
        # Non-aromatic patterns
        non_aromatic_patterns = [
            # "[nH0]1[cH0][nH0][cH0][cH0][cH0]1",
            # "[nH0]1[cH0][nH0][cH0][nH0][cH0]1",
            # "[nH0X2]1[c;R2][c;R2][nH0X3][c;R2][c;R2]1",
            "n1[cH0;$([cH0]=O)]nccc1",
            "n1ccnc[cH0;$([cH0]=O)]1",
            "n1ccn[cH0;$([cH0]=O)][cH0;$([cH0]=O)]1",
            "[nH0]1[cH0;$([cH0]=O)][nH0]cn[cH0;$([cH0]=O)]1",
            "c1ccco1",
            "n1[nH0][cH0;$([cH0]=O)]ccc1",
            "[cH0]1[cH0]cc[cH0]cc1",
            "[nH0;R2]1ccccc1",
            "n1[cH0;$([cH0]=*)]cccc1",
            "n1c[cH0;$([cH0]=*)]ccc1",
            "n1cc[cH0;$([cH0]=*)]cc1",
            "n1ccc[cH0;$([cH0]=*)]c1",
            "n1cccc[cH0;$([cH0]=*)]1",
            "n1cn[cH0;$([cH0]=O)]cc1",
        ]

        # Aromatic patterns
        aromatic_patterns = [
            "c1ccccc1",
            "c1ncncn1",
            "n1ccccc1",
            "n1ccncc1",
            "[nH0;R1]1[c;R1][nH0]ccc1",
            "[nH0]1[nH0]cccc1",
            "n1ccccc1",
            "[nX2]1[nX2]c[nX2]cc1",
            "[nH0]1c[nH0]ccc1",
        ]

        # Ring as set to compare
        ring_set = set(ring)

        # Check if the ring matches any non aromatic pattern
        for smarts in non_aromatic_patterns:
            pattern = Chem.MolFromSmarts(smarts)

            for match in mol.GetSubstructMatches(pattern):
                match_set = set(match)
                if match_set == ring_set:
                    return False

        # Check if the ring matches any aromatic pattern
        for smarts in aromatic_patterns:
            pattern = Chem.MolFromSmarts(smarts)

            for match in mol.GetSubstructMatches(pattern):
                match_set = set(match)
                if match_set == ring_set:
                    return True

        # Special case: all are carbon atoms
        if all(mol.GetAtomWithIdx(idx).GetSymbol() == "C" for idx in ring):
            return True

        # Default case, we consider the ring as non aromatic
        return False
