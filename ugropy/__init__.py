"""ugropy library.

A Python library designed to swiftly and effortlessly obtain the UNIFAC-like
groups from molecules by their names and subsequently integrate them into
inputs for thermodynamic libraries. UNIFAC, PSRK, Dortmund-UNIFAC,
Abdulelah-Gani, and Joback models are implemented.
"""

from .constants import ureg
from .core import instantiate_mol_object
from .core.frag_classes.abdulelah_gani.abdulelah_gani_pst import (
    AbdulelahGaniPSTModel,
)
from .core.frag_classes.abdulelah_gani.abdulelah_gani_pst_result import (
    AGaniPSTFragmentationResult,
)
from .core.frag_classes.base.fragmentation_model import (
    FragmentationModel,
)
from .core.frag_classes.base.fragmentation_result import (
    FragmentationResult,
)
from .core.frag_classes.gibbs_model.gibbs_model import GibbsModel
from .core.frag_classes.gibbs_model.gibbs_result import (
    GibbsFragmentationResult,
)
from .core.frag_classes.joback.joback_model import JobackModel
from .core.frag_classes.joback.joback_result import JobackFragmentationResult
from .core.ilp_solvers.default_solver import DefaultSolver
from .core.ilp_solvers.ilp_solver import ILPSolver
from .groups import Groups
from .models.abdulelah_gani_mod import abdulelah_gani
from .models.abdulelah_gani_pmod import abdulelah_gani_p
from .models.abdulelah_gani_smod import abdulelah_gani_s
from .models.abdulelah_gani_tmod import abdulelah_gani_t
from .models.dortmundmod import dortmund
from .models.jobackmod import joback
from .models.psrkmod import psrk
from .models.unifacmod import unifac


__all__ = [
    # =========================================================================
    # Misc
    # -------------------------------------------------------------------------
    "ureg",
    "writers",
    "instantiate_mol_object",
    "Groups",
    # =========================================================================
    # Fragmentation models classes
    # -------------------------------------------------------------------------
    "FragmentationModel",  # Base class
    "AbdulelahGaniPSTModel",
    "GibbsModel",
    "JobackModel",
    # =========================================================================
    # Fragmentation results classes
    # -------------------------------------------------------------------------
    "FragmentationResult",  # Base class
    "AGaniPSTFragmentationResult",
    "GibbsFragmentationResult",
    "JobackFragmentationResult",
    # =========================================================================
    # ILP Solvers
    # -------------------------------------------------------------------------
    "ILPSolver",  # Base class
    "DefaultSolver",
    # =========================================================================
    # Models instances
    # -------------------------------------------------------------------------
    "abdulelah_gani",
    "abdulelah_gani_p",
    "abdulelah_gani_s",
    "abdulelah_gani_t",
    "dortmund",
    "joback",
    "unifac",
    "psrk",
]
