"""FragmentationModel module.

All ugropy models (joback, unifac, psrk) are instances of the
FragmentationModule class.
"""

from abc import ABC, abstractmethod

import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import numpy as np


class FragmentationModel(ABC):
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

    def detect_groups(self, molecule: Chem.rdchem.Mol) -> pd.DataFrame:
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

    @abstractmethod
    def set_fragmentation_result(
        self,
        molecule: Chem.rdchem.Mol,
        subgroups_occurrences: dict,
        subgroups_atoms_indexes: dict,
    ) -> "FragmentationResult":

        raise NotImplementedError("Abstract Method not implemented.")


class FragmentationResult:
    def __init__(self, molecule: Chem.rdchem.Mol, subgroups: dict):
        self.mol_object = molecule
        self.subgroups = subgroups

    def draw(
        mol_object: Chem.rdchem.Mol,
        subgroups: dict,
        model,  # El tipo de model depende de tu implementación
        title: str = "",
        width: int = 400,
        height: int = 200,
        title_font_size: int = 12,
        legend_font_size: int = 12,
        font: str = "Helvetica",
    ) -> Image.Image:
        """Create a PIL image of the fragmentation result with a legend."""

        # Ajustar los subgrupos a los átomos de la molécula
        fit = fit_atoms(mol_object, subgroups, model)

        # Generar los colores para cada subgrupo
        how_many_subgroups = len(fit.keys())
        colors_rgb = _generate_distinct_colors(how_many_subgroups)

        highlight = []
        atoms_colors = {}

        for idx, (subgroup, atoms) in enumerate(fit.items()):
            atms = np.array(atoms).flatten()
            highlight.extend(atms.tolist())

            for at in atms:
                atoms_colors[int(at)] = colors_rgb[idx]

        # Crear la imagen de la molécula usando MolToImage
        img = Draw.MolToImage(
            mol_object,
            size=(width, height),
            highlightAtoms=highlight,
            highlightAtomColors={i: tuple(c[:3]) for i, c in atoms_colors.items()}  # RGB sin el canal alpha
        )

        # Crear una imagen más grande para añadir la leyenda y el título
        new_height = height + (legend_font_size + 5) * how_many_subgroups + 50  # Espacio para la leyenda y título
        final_img = Image.new("RGB", (width, new_height), "white")
        
        # Pegar la imagen de la molécula en la parte superior
        final_img.paste(img, (0, 0))

        # Crear un objeto ImageDraw para dibujar la leyenda y el título
        draw = ImageDraw.Draw(final_img)

        # Cargar la fuente o usar la predeterminada si no está disponible
        try:
            font_title = ImageFont.truetype("arial.ttf", title_font_size)
            font_legend = ImageFont.truetype("arial.ttf", legend_font_size)
        except IOError:
            font_title = ImageFont.load_default()
            font_legend = ImageFont.load_default()

        # Dibujar el título
        draw.text((width / 2, height + 10), title, fill="black", font=font_title, anchor="ms")

        # Dibujar la leyenda
        for i, (subgroup, color) in enumerate(atoms_colors.items()):
            r, g, b = [int(255 * c) for c in color[:3]]  # Convertir a RGB 0-255
            rect_color = (r, g, b)

            # Dibujar un cuadro de color para cada subgrupo
            draw.rectangle([10, height + 50 + i * (legend_font_size + 5), 30, height + 50 + (i + 1) * (legend_font_size + 5)], fill=rect_color)
            
            # Dibujar el nombre del subgrupo
            draw.text((40, height + 50 + i * (legend_font_size + 5)), f"Subgrupo {i + 1}", fill="black", font=font_legend)

        return final_img

# Generar colores distintos (como antes)
def _generate_distinct_colors(n: int) -> list:
    base_colors = np.array(
        [
            [0.12156863, 0.46666667, 0.70588235],  # azul
            [1.0, 0.49803922, 0.05490196],  # naranja
            [0.17254902, 0.62745098, 0.17254902],  # verde
            [0.83921569, 0.15294118, 0.15686275],  # rojo
            [0.58039216, 0.40392157, 0.74117647],  # púrpura
            [0.54901961, 0.3372549, 0.29411765],  # marrón
            [0.89019608, 0.46666667, 0.76078431],  # rosa
            [0.49803922, 0.49803922, 0.49803922],  # gris
            [0.7372549, 0.74117647, 0.13333333],  # amarillo
            [0.09019608, 0.74509804, 0.81176471],  # cian
        ]
    )
    colors = [base_colors[i % len(base_colors)] for i in range(n)]
    return [(color[0], color[1], color[2], 0.65) for color in colors]

# Función de prueba
mol = Chem.MolFromSmiles('CCO')
subgroups = {'Grupo 1': [0, 1], 'Grupo 2': [2]}  # Ejemplo
img_with_legend = draw(mol, subgroups, model=None, title="Molécula con subgrupos")
img_with_legend.show()
