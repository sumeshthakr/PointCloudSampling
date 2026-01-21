"""Poisson disk sampling (3D) - refactored and optimized.

Public API: Poisson
"""
from __future__ import annotations

from typing import List, Tuple
import numpy as np


class Poisson:
    """Poisson disk sampling in a 3D rectangular box.

    Attributes:
        k: int - attempts per active point
        r: float - minimum allowed distance between points
        max_x, max_y, max_z: float - box extents
    """

    def __init__(self, k: int, r: float, max_x: float, max_y: float, max_z: float):
        self.k = int(k)
        self.r = float(r)
        self.max_x = float(max_x)
        self.max_y = float(max_y)
        self.max_z = float(max_z)

        # cell size: a smaller cell than r ensures at most one point per cell
        self.a = self.r / np.sqrt(3)

        # number of voxels along each axis
        self.nx = int(self.max_x / self.a) + 1
        self.ny = int(self.max_y / self.a) + 1
        self.nz = int(self.max_z / self.a) + 1

        self.reset_voxels()

    def reset_voxels(self) -> None:
        """Initialize voxel grid and sample list."""
        self.voxels = {}
        self.samples: List[np.ndarray] = []

    def get_cell_coords(self, point: Tuple[float, float, float]) -> Tuple[int, int, int]:
        return int(point[0] // self.a), int(point[1] // self.a), int(point[2] // self.a)

    def get_neighbour(self, coords: Tuple[int, int, int]) -> List[int]:
        """Return indices of samples in neighbouring voxels (within +/-2 in x,y and +/-1 in z).

        This covers the region that may contain points closer than r.
        """
        neighbours: List[int] = []
        cx, cy, cz = coords
        for dx in range(-2, 3):
            nx = cx + dx
            if nx < 0 or nx >= self.nx:
                continue
            for dy in range(-2, 3):
                ny = cy + dy
                if ny < 0 or ny >= self.ny:
                    continue
                for dz in range(-2, 3):
                    nz = cz + dz
                    if nz < 0 or nz >= self.nz:
                        continue
                    idx = self.voxels.get((nx, ny, nz))
                    if idx is not None:
                        neighbours.append(idx)
        return neighbours

    def point_valid(self, pt: Tuple[float, float, float]) -> bool:
        """Check whether pt is at least r away from existing samples (vectorized)."""
        cell_coords = self.get_cell_coords(pt)
        ids = self.get_neighbour(cell_coords)
        if not ids:
            return True
        nearby = np.vstack([self.samples[i] for i in ids])
        diff = nearby - np.asarray(pt)
        d2 = np.sum(diff * diff, axis=1)
        return not np.any(d2 < (self.r * self.r))

    def get_point(self, refpt: Tuple[float, float, float]):
        """Try up to k times to generate a valid point in the annulus [r, 2r] around refpt."""
        for _ in range(self.k):
            rho = np.random.uniform(self.r, 2 * self.r)
            theta = np.random.uniform(0, 2 * np.pi)
            phi = np.random.uniform(0, np.pi)
            x = refpt[0] + rho * np.sin(phi) * np.cos(theta)
            y = refpt[1] + rho * np.sin(phi) * np.sin(theta)
            z = refpt[2] + rho * np.cos(phi)
            if not (0 <= x < self.max_x and 0 <= y < self.max_y and 0 <= z < self.max_z):
                continue
            candidate = (x, y, z)
            if self.point_valid(candidate):
                return candidate
        return None

    def sample(self) -> np.ndarray:
        """Run Poisson sampling algorithm and return points as an (N,3) numpy array."""
        self.reset_voxels()
        pt = (np.random.uniform(0, self.max_x), np.random.uniform(0, self.max_y), np.random.uniform(0, self.max_z))
        self.samples = [np.array(pt)]
        self.voxels[self.get_cell_coords(pt)] = 0
        active = [0]

        while active:
            idx = int(np.random.choice(active))
            refpt = tuple(self.samples[idx])
            pt = self.get_point(refpt)
            if pt is not None:
                self.samples.append(np.array(pt))
                ns = len(self.samples) - 1
                active.append(ns)
                self.voxels[self.get_cell_coords(pt)] = ns
            else:
                active.remove(idx)

        return np.vstack(self.samples)
