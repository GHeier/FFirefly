"""
Many-body physics utilities for FFirefly.

This module provides tools for many-body calculations including:
- Interface with TRIQS library (interface_triqs)
- Loading TRIQS Hamiltonians (load_triqs_H)
"""

from . import interface_triqs
from . import load_triqs_H

# Export commonly used functions for convenience
from .interface_triqs import (
    fill_triqs_from_field,
)

from .load_triqs_H import get_energy_mesh

__all__ = [
    'interface_triqs',
    'load_triqs_H',
    'fill_triqs_from_field',
    'get_energy_mesh'
]
