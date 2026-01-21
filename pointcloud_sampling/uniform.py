"""Uniform sampling utilities."""
from __future__ import annotations

from typing import Tuple
import numpy as np


def sample_uniform(n: int, x_range: Tuple[float, float], y_range: Tuple[float, float], z_range: Tuple[float, float]) -> np.ndarray:
    """Return n points sampled uniformly in the given box ranges as (n,3) array."""
    xs = np.random.uniform(x_range[0], x_range[1], size=n)
    ys = np.random.uniform(y_range[0], y_range[1], size=n)
    zs = np.random.uniform(z_range[0], z_range[1], size=n)
    return np.column_stack([xs, ys, zs])
