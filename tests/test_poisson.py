import pytest
import numpy as np
from src.point_cloud_sampling.poisson import PoissonSampler

def test_poisson_bounds():
    max_dims = (10, 10, 10)
    sampler = PoissonSampler(max_attempts=30, radius=2.0, max_x=max_dims[0], max_y=max_dims[1], max_z=max_dims[2])
    points = sampler.sample()
    points_np = np.array(points)
    
    assert len( points ) > 0
    assert np.all(points_np >= 0)
    assert np.all(points_np < max_dims)

def test_poisson_min_distance():
    radius = 1.5
    sampler = PoissonSampler(max_attempts=30, radius=radius, max_x=10, max_y=10, max_z=10)
    points = sampler.sample()
    points_np = np.array(points)
    
    # Calculate pairwise distances
    diff = points_np[:, np.newaxis, :] - points_np[np.newaxis, :, :]
    dist_sq = np.sum(diff**2, axis=-1)
    
    # Fill diagonal with infinity to ignore self-distance
    np.fill_diagonal(dist_sq, np.inf)
    
    min_dist_sq = np.min(dist_sq)
    assert min_dist_sq >= radius**2 - 1e-9 # allow small float error

def test_poisson_empty():
    # Should handle small bounds/large radius gracefully
    sampler = PoissonSampler(30, 20, 5, 5, 5) # Radius > bounds
    points = sampler.sample()
    # Should assume at least one point fits if bounds > 0? 
    # Radius 20 in 5x5x5 box. First point fits. neighbors will be invalid.
    assert len(points) == 1
