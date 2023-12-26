"""Writers module.

The witers module contains functions to construct the neccesary inputs for
other thermodynamics libraries.

Supported:

- Clapeyron.jl: https://github.com/ClapeyronThermo/Clapeyron.jl
"""
from .clapeyron import to_clapeyron


__all__ = ["to_clapeyron"]
