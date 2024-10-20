"""FragmentationResult class.

Base class to set fragmentation results. implements the drawing methods of the
solutions.
"""

import importlib

import numpy as np

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


class FragmentationResult:
    """Fragmentation result class.

    This class is used to store the fragmentation results and provide methods
    to visualize them.

    Parameters
    ----------
    molecule : Chem.rdchem.Mol
        Molecule to fragment.
    subgroups : dict
        Dictionary with the subgroups and the number of times they appear in
        the molecule.
    subgroups_atoms_indexes : dict
        Dictionary with the subgroups and the atoms indexes that belong to
        each subgroup.

    Attributes
    ----------
    molecule : Chem.rdchem.Mol
        Molecule to fragment.
    subgroups : dict
        Dictionary with the subgroups and the number of times they appear in
        the molecule.
    subgroups_atoms : dict
        Dictionary with the subgroups and the atoms indexes that belong to
        each subgroup.
    """

    def __init__(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups: dict,
        subgroups_atoms_indexes: dict,
    ):
        self.molecule = molecule
        self.subgroups = subgroups
        self.subgroups_atoms = subgroups_atoms_indexes

    def get_solution_svg(
        self,
        title: str = "",
        width: float = 400,
        height: float = 200,
        title_font_size: float = 12,
        legend_font_size: float = 12,
        font: str = "Helvetica",
    ) -> str:
        """Create a SVG figure string of the fragmentation result.

        Parameters
        ----------
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
            self.molecule,
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

        for i, (name, color) in enumerate(subgroups_colors.items()):
            r, g, b, _ = color
            rect_color = f"rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})"

            legend_name = name.replace("<", "&lt;").replace(">", "&gt;")

            # Color rectangle for the legend of the group
            legend += (
                f'<rect x="1" y="{5 + i * 25}" '
                f'width="{legend_font_size * 1.5}" '
                f'height="{legend_font_size * 1.5}" fill="{rect_color}" />'
            )

            # Group name and occurences for the legend
            legend += (
                f'<text x="{legend_font_size * 1.6}" y="{20 + i * 25}" '
                f'font-family="{font}" font-size="{legend_font_size}" '
                f'fill="black">{legend_name}: {self.subgroups[name]}</text>'
            )

        # =====================================================================
        # Create title
        # =====================================================================
        title = (
            f'<text x="{width/2}" y="40" '
            f'font-family="{font}" font-size="{title_font_size}" '
            f'font-weight="bold" fill="black" '
            f'text-anchor="middle">{title}</text>'
        )

        # =====================================================================
        # Set title and legend to figure
        # =====================================================================
        svg_complete = svg.replace("</svg>", f"{legend}{title}</svg>")

        return svg_complete

    def draw(
        self,
        title: str = "",
        width: float = 400,
        height: float = 200,
        title_font_size: float = 12,
        legend_font_size: float = 12,
        font: str = "Helvetica",
    ):
        """Create a IPython SVG object of the fragmentation result.

        This function is meant to be used in Jupyter notebooks to directly
        obtain the visual SVG figure on the notebook. It requires the IPython
        library to be installed.

        Parameters
        ----------
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
        IPython.display.SVG
            SVG object of the fragmentation result.

        Raises
        ------
        ImportError
            IPython is not installed.
        """
        svg = self.get_solution_svg(
            title, width, height, title_font_size, legend_font_size, font
        )

        try:
            svg_lib = importlib.import_module("IPython.display")
        except ImportError:
            raise ImportError(
                "IPython is required to display the SVG. This function is"
                "meant to be used in Jupyter notebooks."
            )

        return svg_lib.SVG(svg)

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
