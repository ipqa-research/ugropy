"""Fragmentation module."""

from typing import List, Union

from rdkit import Chem

from ugropy.core.draw_molecule import draw
from ugropy.fragmentation_models.fragmentation_model import FragmentationModel


class Fragmentation:
    """Fragmentation object. Is the return of the get_groups function.

    Parameters
    ----------
    mol_subgroups : Union[dict, List[dict]]
        Subgroups of mol_object.
    mol_object : Chem.rdchem.Mol
        RDKit Mol object.
    model: FragmentationModel
        FragmentationModel object.

    Attributes
    ----------
    subgroups : dict
        Subgroups of the molecule.
    mol_object : Chem.rdchem.Mol
        RDKit Mol object.
    """

    def __init__(
        self,
        mol_subgroups: Union[dict, List[dict]],
        mol_object: Chem.rdchem.Mol,
        model: FragmentationModel,
    ) -> None:
        self.subgroups = mol_subgroups
        self.mol_object = mol_object
        self._model = model

    def draw(
        self,
        title: str = "",
        width: float = 400,
        height: float = 200,
        title_font_size: float = 12,
        legend_font_size: float = 12,
        font: str = "Helvetica",
    ) -> Union[str, List[str]]:
        """Create a svg representation of the fragmentation result.

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
        Union[str, List[str]]
            SVG of the fragmentation solution/s.
        """
        if not isinstance(self.subgroups, dict):
            svg_list = []
            for sol in self.subgroups:
                svg = draw(
                    self.mol_object,
                    sol,
                    self._model,
                    title,
                    width,
                    height,
                    title_font_size,
                    legend_font_size,
                    font,
                )

                svg_list.append(svg)
            return svg_list
        else:
            svg = draw(
                self.mol_object,
                self.subgroups,
                self._model,
                title,
                width,
                height,
                title_font_size,
                legend_font_size,
                font,
            )
            return svg
