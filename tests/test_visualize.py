import pytest
import numpy as np
import matplotlib
matplotlib.use('Agg')

from src.point_cloud_sampling.visualize import plot_samples, plot_comparison, plot_density_comparison


@pytest.fixture
def sample_points():
    np.random.seed(0)
    return np.random.uniform(0, 10, (50, 3))


@pytest.fixture
def two_sample_sets():
    np.random.seed(0)
    poisson = np.random.uniform(0, 10, (50, 3))
    uniform = np.random.uniform(0, 10, (50, 3))
    return poisson, uniform


def test_plot_samples_save(tmp_path, sample_points):
    path = str(tmp_path / "test_plot.png")
    plot_samples(sample_points, title="Test", save_path=path)
    import os
    assert os.path.exists(path)
    assert os.path.getsize(path) > 0


def test_plot_samples_empty():
    empty = np.array([]).reshape(0, 3)
    # Should not raise, just print a message
    plot_samples(empty)


def test_plot_comparison_save(tmp_path, two_sample_sets):
    poisson, uniform = two_sample_sets
    path = str(tmp_path / "test_comparison.png")
    plot_comparison(poisson, uniform, save_path=path)
    import os
    assert os.path.exists(path)
    assert os.path.getsize(path) > 0


def test_plot_density_comparison_save(tmp_path, two_sample_sets):
    poisson, uniform = two_sample_sets
    path = str(tmp_path / "test_density.png")
    plot_density_comparison(poisson, uniform, save_path=path)
    import os
    assert os.path.exists(path)
    assert os.path.getsize(path) > 0
