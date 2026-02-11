# Point Cloud Sampling

<p align="center">
  <img src="logo.svg" alt="Point Cloud Sampling Logo" width="600"/>
</p>


A high-performance Python package for 3D point cloud sampling, implementing **Poisson Disk Sampling** and **Uniform Sampling**.

**Now optimized with Numba!** 🚀

## Features

- **Poisson Disk Sampling**: Generates points with a minimum distance constraint ($r$) using Bridson's algorithm.
    - **Optimized**: Uses `numba` for JIT compilation, achieving **~140x speedup** over pure Python.
- **Uniform Sampling**: Generates points uniformly distributed within a bounding box.
- **Visualization**: Built-in 3D visualization using Matplotlib.

## Performance

Benchmark results for generating ~22,000 points in a $50 \times 50 \times 50$ volume with $r=1.5$:

| Implementation | Time (avg) | Points/sec |
| :--- | :--- | :--- |
| Pure Python (Previous) | ~27.0 s | ~800 |
| **Numba Optimized** | **0.19 s** | **~118,000** |

## Installation

This package requires Python 3.8+ and `numba`.

### From PyPI (Coming Soon)

```bash
pip install point-cloud-sampling
```

### From Source

```bash
# Clone the repository
git clone https://github.com/sumeshthakr/PointCloudSampling.git
cd PointCloudSampling

# Create a virtual environment (recommended)
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies (includes numba)
pip install .
```

## Usage

### Poisson Sampling

Poisson Disk Sampling is ideal for generating blue noise sample patterns.

```python
from point_cloud_sampling import PoissonSampler, plot_samples

# Initialize sampler
# max_attempts: Limit of attempts to find a neighbor (standard is 30)
# radius: Minimum distance between points
# max_x, max_y, max_z: Bounding box dimensions
sampler = PoissonSampler(max_attempts=30, radius=1.5, max_x=10, max_y=10, max_z=10)

# Generate samples (First run may take a moment to compile)
points = sampler.sample()

print(f"Generated {len(points)} points.")

# Visualize
# plot_samples(points)
```

### Uniform Sampling

```python
from point_cloud_sampling import UniformSampler, plot_samples

# Generate 100 points
points = UniformSampler.sample(
    n_points=100,
    x_range=(0, 10),
    y_range=(0, 10),
    z_range=(0, 10)
)

plot_samples(points)
```

## Future Plans

-   **Parallelization**: Implement parallel processing for even faster generation on multi-core CPUs.
-   **GPU Support**: Explore CUDA support using `numba.cuda` for massive point cloud generation.
-   **More Shapes**: Support sampling within spheres, cylinders, and custom meshes.
-   **Blue Noise Properties**: Add tools to analyze the spectral properties of the generated samples.
-   **Export Formats**: Support exporting to `.ply`, `.xyz`, and `.pcd` formats.

## License

MIT License. See [LICENSE](LICENSE) for details.
