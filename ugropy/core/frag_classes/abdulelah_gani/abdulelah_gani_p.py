"""Joback fragmentation module."""

from typing import List, Union

import pandas as pd

from rdkit import Chem

from ugropy.core.frag_classes.abdulelah_gani.abdulelah_gani_p_result import (
    AGaniPFragmentationResult,
)
from ugropy.core.frag_classes.base.fragmentation_model import (
    FragmentationModel,
)
from ugropy.core.ilp_solvers.default_solver import DefaultSolver
from ugropy.core.ilp_solvers.ilp_solver import ILPSolver


class AbdulelahGaniPrimaryModel(FragmentationModel):
    """Abdulelah-Gani model dedicated to properties estimation models.

    Class to construct the primary structures detector for the Abdulelah-Gani
    properties estimation model :cite:p:`gani`.

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
    ) -> None:

        super().__init__(
            subgroups=subgroups,
            allow_overlapping=False,
            fragmentation_result=AGaniPFragmentationResult,
        )
        self.subgroups_info = subgroups_info

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        solver: ILPSolver = DefaultSolver,
        search_multiple_solutions: bool = False,
    ) -> Union[AGaniPFragmentationResult, List[AGaniPFragmentationResult]]:
        """Get the groups of a molecule.

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
            subgroups_info=self.subgroups_info,
        )

        return sol

    def mol_preprocess(self, mol: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
        """

        Parameters
        ----------
        mol : Chem.rdchem.Mol
            _description_

        Returns
        -------
        Chem.rdchem.Mol
            _description_
        """

        # Clone of the molecule to kekulize
        kekulized_mol = Chem.Mol(mol)
        Chem.Kekulize(kekulized_mol, clearAromaticFlags=True)

        # Step 1: Identify all ring on the molecule
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()

        # Step 2: Get only the aromatic rings
        aromatic_rings = [
            ring
            for ring in atom_rings
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]
        
        true_aromatic_rings = []
        
        # Step 3: Preprocess all rings that are not exclusively carbon atoms
        #         to make them aliphatic.
        for ring in aromatic_rings:
            # All ring's atoms are carbon?
            if not all(
                mol.GetAtomWithIdx(idx).GetSymbol() == "C" for idx in ring
            ):  
                # Specific rings adhocs that should be kept aromatic  
                if is_aromatic(mol, ring):
                    true_aromatic_rings.append(ring)
                    continue
                
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
            else:
                true_aromatic_rings.append(ring)

        # Step 4: Revisit all the Carbon-only aromatic rings and set them as
        #         aromatic. This is done because the previous step may have
        #         marked some atoms as non-aromatic when two rings share atoms.
        for ring in true_aromatic_rings:
            # Restore atoms aromaticity
            for idx in ring:
                mol.GetAtomWithIdx(idx).SetIsAromatic(True)

            # Restore bonds aromaticity
            for i in range(len(ring)):
                atom1, atom2 = ring[i], ring[(i + 1) % len(ring)]
                bond = mol.GetBondBetweenAtoms(atom1, atom2)
                bond.SetBondType(Chem.rdchem.BondType.AROMATIC)

        return mol

def is_aromatic(mol, ring):
    """
    Verifica si un anillo coincide con un patrón SMARTS en la molécula original.

    Parameters
    ----------
    mol : Chem.Mol
        Molécula original.
    ring : list[int]
        Lista de índices de los átomos del anillo.

    Returns
    -------
    bool
        True si el anillo coincide con algún patrón SMARTS, False en caso contrario.
    """
    non_aromatic_patterns = [
        "[nH0]1[cH0][nH0][cH0][cH0][cH0]1"
    ]
    
    # Lista de patrones SMARTS
    aromatic_patterns = [
        "c1ncncn1",
        "n1ccccc1",
        "n1ccncc1",
        "[nH0]1c[nH0]ccc1",
        "[nH0]1[nH0]cccc1",
        "n1ccccc1"
    ]

    # Convertir la lista de índices del anillo en un set para comparación eficiente
    ring_set = set(ring)
    
    for smarts in non_aromatic_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            raise ValueError(f"El patrón SMARTS '{smarts}' no es válido.")
        
        # Buscar coincidencias del patrón en la molécula
        for match in mol.GetSubstructMatches(pattern):
            # Convertir la coincidencia a set para comparación
            match_set = set(match)
            if match_set == ring_set:
                return False  # El patrón coincide exactamente con el anillo

    for smarts in aromatic_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            raise ValueError(f"El patrón SMARTS '{smarts}' no es válido.")
        
        # Buscar coincidencias del patrón en la molécula
        for match in mol.GetSubstructMatches(pattern):
            # Convertir la coincidencia a set para comparación
            match_set = set(match)
            if match_set == ring_set:
                return True  # El patrón coincide exactamente con el anillo
        
    return False  # Ningún patrón coincide
