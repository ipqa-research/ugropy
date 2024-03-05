"""draw_molecule module."""

from typing import List

import numpy as np

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from ugropy.core.fit_atoms_indexes import fit_atoms
from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


def draw(
    mol_object: Chem.rdchem.Mol,
    subgroups: dict,
    model: FragmentationModel,
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
    # Fit subgroups into the mol's atoms
    # =====================================================================
    fit = fit_atoms(mol_object, subgroups, model)

    # =====================================================================
    # Generate the colors for each subgroup (max 10 TODO)
    # =====================================================================
    how_many_subgroups = len(fit.keys())

    colors_rgb = _generate_distinct_colors(how_many_subgroups)

    highlight = []
    atoms_colors = {}
    subgroups_colors = {}

    for idx, (subgroup, atoms) in enumerate(fit.items()):
        atms = np.array(atoms).flatten()

        subgroups_colors[subgroup] = colors_rgb[idx]
        highlight.extend(atms.tolist())

        for at in atms:
            atoms_colors[int(at)] = colors_rgb[idx]

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(
        mol_object,
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


def _generate_distinct_colors(n: int) -> List[tuple]:
    """Return distinguishable colors in (r,g,b,a) format (10 different max).

    Parameters
    ----------
    n : int
        Number of colors desired.

    Returns
    -------
    List[tuple]
        Colors.
    """
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
            [0.09019608, 0.74509804, 0.81176471],  # cyan
        ]
    )
    colors = [base_colors[i % len(base_colors)] for i in range(n)]
    colors_rgb = [(color[0], color[1], color[2], 0.65) for color in colors]
    return colors_rgb
