import pytest
import numpy as np
from src.point_cloud_sampling.uniform import UniformSampler

def test_uniform_bounds():
    n = 100
    x_range = (0, 10)
    y_range = (5, 15)
    z_range = (-5, 5)
    
    points = UniformSampler.sample(n, x_range, y_range, z_range)
    
    assert len(points) == n
    assert np.all(points[:, 0] >= x_range[0])
    assert np.all(points[:, 0] < x_range[1])
    assert np.all(points[:, 1] >= y_range[0])
    assert np.all(points[:, 1] < y_range[1])
    assert np.all(points[:, 2] >= z_range[0])
    assert np.all(points[:, 2] < z_range[1])
