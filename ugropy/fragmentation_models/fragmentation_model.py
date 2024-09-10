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

    def __init__(self, subgroups: pd.DataFrame) -> None:
        self.subgroups = subgroups

        # Instantiate all de mol object from their smarts representation
        detection_mols = {}

        for group, row in self.subgroups.iterrows():
            detection_mols[group] = Chem.MolFromSmarts(row["smarts"])

    def get_groups(
        self,
        identifier: Union[str, Chem.rdchem.Mol],
        identifier_type: str = "name",
        ilp_solver: str = "cbc",
    ) -> "FragmentationResult":
        
        # RDKit Mol object
        mol_object = instantiate_mol_object(identifier, identifier_type)

        # Direct detection of groups presence and occurences
        detections = self._detect_groups(mol_object)

        # First return
        if detections == {}:  # No groups detected
            return self.set_fragmentation_result(mol_object, {}, {})
        
        # Check overlapping groups
        has_overlap, overlapping_atoms = check_has_overlapping_groups(
            mol_object, detections
        )

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


    def _detect_groups(self, molecule: Chem.rdchem.Mol) -> dict:
        """Detect all the groups in the molecule.

        Return a dictionary with the detected groups as keys and a tuple of
        tuples containing the molecule's atoms that participate in the group
        occurrences.

        Parameters
        ----------
        mol : Chem.rdchem.Mol
            Molecule to detect the groups.

        Returns
        -------
        dict
            Detected groups in the molecule.
        """
        detected_groups = {}
        for group, mol in self.detection_mols.items():
            matches = molecule.GetSubstructMatches(mol)

            if matches:
                detected_groups[group] = matches

        return detected_groups


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
