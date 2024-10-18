import sys
import os

from IPython.display import display, HTML

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from ugropy.core.get_rdkit_object import instantiate_mol_object


# path tests
project_root = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "../../../tests/")
)

sys.path.append(project_root)


# Display function
def display_case_module(cases: list):
    for case in cases:
        mol = instantiate_mol_object(case.identifier, case.identifier_type)

        drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        svg = drawer.GetDrawingText().replace("svg:", "")

        # HTML block
        case_html = f"""
        <div style="border: 1px solid #ccc; padding: 10px; margin-bottom: 20px; border-radius: 5px;">
            <h3 style="color: #2C3E50;">Identifier: {case.identifier}</h3>
            <p><strong style="color: #3498DB;">UNIFAC Result:</strong> {case.unifac_result}</p>
            <p><strong style="color: #3498DB;">PSRK Result:</strong> {case.psrk_result}</p>
            <p><strong style="color: #3498DB;">Joback Result:</strong> {case.joback_result}</p>
            <div style="text-align: center;">
                {svg}
            </div>
        </div>
        """

        # display HTML
        display(HTML(case_html))
