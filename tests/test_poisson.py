import numpy as np
from pointcloud_sampling import Poisson


def test_poisson_min_distance():
    p = Poisson(k=30, r=1.0, max_x=10, max_y=10, max_z=10)
    pts = p.sample()
    assert pts.shape[1] == 3
    if pts.shape[0] > 1:
        diffs = pts[:, None, :] - pts[None, :, :]
        d2 = np.sum(diffs * diffs, axis=-1)
        np.fill_diagonal(d2, np.inf)
        min_d2 = d2.min()
        assert min_d2 >= (0.99 * p.r) ** 2
