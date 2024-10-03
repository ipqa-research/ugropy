"""FragmentationModel module.

All ugropy models (joback, unifac, psrk) are instances of the
FragmentationModule class.
"""

from typing import Union

import pandas as pd

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

import numpy as np

from ugropy.core.checks import check_has_overlapping_groups
from ugropy.core.get_rdkit_object import instantiate_mol_object


class FragmentationModel:
    """FragmentationModel class.

    All ugropy supported models are an instance of this class. This class must
    be inherited to create a new type of FragmentationModel.

    Parameters
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct).

    Attributes
    ----------
    subgroups : pd.DataFrame
        Model's subgroups. Index: 'group' (subgroups names). Mandatory columns:
        'smarts' (SMARTS representations of the group to detect its precense in
        the molecule), 'molecular_weight' (molecular weight of the subgroups
        used to check that the result is correct).
    detection_mols : dict
        Dictionary cotaining all the rdkit Mol object from the detection_smarts
        subgroups column.
    """

    def __init__(
        self,
        subgroups: pd.DataFrame,
        allow_overlapping: bool = False,
        check_molecular_weight: bool = False,
    ) -> None:
        self.subgroups = subgroups
        self.allow_overlapping = allow_overlapping
        self.check_molecular_weight = check_molecular_weight

        # Instantiate all de mol object from their smarts representation
        self.detection_mols = {}

        for group, row in self.subgroups.iterrows():
            self.detection_mols[group] = (Chem.MolFromSmarts(row["smarts"]))

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        ilp_solver: str = "cbc",
    ) -> "FragmentationResult":

        # RDKit Mol object
        mol_object = instantiate_mol_object(identifier, identifier_type)

        # Direct detection of fragments presence and its atoms indexes
        detections = self.detect_fragments(mol_object)

        # First return
        if detections == {}:  # No groups detected
            return self.set_fragmentation_result(mol_object, {}, {})

        # Check overlapping groups
        has_overlap, overlapping_atoms = check_has_overlapping_groups(
            mol_object, detections
        )

        # Second return
        if not has_overlap:
            return self.set_fragmentation_result(mol_object, detections, overlapping_atoms)


    def set_fragmentation_result(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups_occurrences: dict,
        subgroups_atoms_indexes: dict,
    ) -> "FragmentationResult":

        result = FragmentationResult(
            molecule, subgroups_occurrences, subgroups_atoms_indexes
        )

        return result

    def detect_fragments(self, molecule: Chem.rdchem.Mol) -> dict:
        """Detect all the fragments in the molecule.

        Return a dictionary with the detected fragments as keys and a tuple
        with the atoms indexes of the fragment as values. For example, n-hexane
        for the UNIFAC model will return:
        
        {
            'CH3_0': (0,),
            'CH3_1': (5,),
            'CH2_0': (1,),
            'CH2_1': (2,),
            'CH2_2': (3,),
            'CH2_3': (4,)
        }

        You may note that multiple occurrence of a fragment name will be 
        indexed. The convention is always: <fragment_name>_i where `i` is the
        index of the occurrence.

        Parameters
        ----------
        mol : Chem.rdchem.Mol
            Molecule to detect the fragments.

        Returns
        -------
        dict
            Detected fragments in the molecule.
        """
        detected_fragments = {}
        
        for fragment_name, mol in self.detection_mols.items():
            matches = molecule.GetSubstructMatches(mol)

            if matches:
                for i, atoms_tuple in enumerate(matches):
                    detected_fragments[f"{fragment_name}_{i}"] = atoms_tuple

        return detected_fragments


class FragmentationResult:
    def __init__(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
    ):
        self.mol_object = molecule
        self.subgroups = subgroups
        self.subgroups_atoms = subgroups_atoms_indexes

    def draw(
        self,
        title: str = "",
        width: float = 400,
        height: float = 200,
        title_font_size: float = 12,
        legend_font_size: float = 12,
        font: str = "Helvetica",
    ) -> str:
        """Create a svg representation of the fragmentation result.

        Parameters
        ----------
        mol_object : Chem.rdchem.Mol
            RDKit Mol object.
        mol_subgroups : Union[dict, List[dict]]
            Subgroups of mol_object.
        model: FragmentationModel
            FragmentationModel object.
        title : str, optional
            Graph title, by default ""
        width : int, optional
            Graph width, by default 400
        height : int, optional
            Graph height, by default 200
        title_font_size : int, optional
            Font size of graph's title, by default 12
        legend_font_size : int, optional
            Legend font size, by default 12
        font : str, optional
            Text font, by default "Helvetica"

        Returns
        -------
        str
            SVG string.
        """
        # =====================================================================
        # Generate the colors for each subgroup (max 10 TODO)
        # =====================================================================
        how_many_subgroups = len(self.subgroups_atoms.keys())

        colors_rgb = self._generate_distinct_colors(how_many_subgroups)

        highlight = []
        atoms_colors = {}
        subgroups_colors = {}

        for idx, (subgroup, atoms) in enumerate(self.subgroups_atoms.items()):
            atms = np.array(atoms).flatten()

            subgroups_colors[subgroup] = colors_rgb[idx]
            highlight.extend(atms.tolist())

            for at in atms:
                atoms_colors[int(at)] = colors_rgb[idx]

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(
            self.mol_object,
            highlightAtoms=highlight,
            highlightBonds=[],
            highlightAtomColors=atoms_colors,
        )
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")

        # =====================================================================
        # Create legend
        # =====================================================================
        legend = ""
        i = 0
        for i, (name, color) in enumerate(subgroups_colors.items()):
            r, g, b, _ = color
            rect_color = f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})"
            name = name.replace("<", "&lt;").replace(">", "&gt;")
            legend += f'<rect x="1" y="{5 + i * 25}" width="{legend_font_size * 1.5}" height="{legend_font_size * 1.5}" fill="{rect_color}" />'  # noqa
            legend += f'<text x="{legend_font_size * 1.6}" y="{20 + i * 25}" font-family="{font}" font-size="{legend_font_size}" fill="black">{name}</text>'  # noqa

        # =====================================================================
        # Create title
        # =====================================================================
        title = f'<text x="{width/2}" y="40" font-family="{font}" font-size="{title_font_size}" font-weight="bold" fill="black" text-anchor="middle">{title}</text>'  # noqa

        # =====================================================================
        # Set title and legend to figure
        # =====================================================================
        svg_with_legend = svg.replace("</svg>", f"{legend}{title}</svg>")

        return svg_with_legend

    def _generate_distinct_colors(self, n: int) -> list:
        base_colors = np.array(
            [
                [0.12156863, 0.46666667, 0.70588235],  # blue
                [1.0, 0.49803922, 0.05490196],  # orange
                [0.17254902, 0.62745098, 0.17254902],  # green
                [0.83921569, 0.15294118, 0.15686275],  # red
                [0.58039216, 0.40392157, 0.74117647],  # purple
                [0.54901961, 0.3372549, 0.29411765],  # brown
                [0.89019608, 0.46666667, 0.76078431],  # pink
                [0.49803922, 0.49803922, 0.49803922],  # gray
                [0.7372549, 0.74117647, 0.13333333],  # yellow
                [0.09019608, 0.74509804, 0.81176471],  # cian
            ]
        )
        colors = [base_colors[i % len(base_colors)] for i in range(n)]
        return [(color[0], color[1], color[2], 0.65) for color in colors]
