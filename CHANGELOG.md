# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - 2026-02-11

### Added
- **Visualization module** (`visualize.py`) with dedicated comparison plotting functions:
  - `plot_samples()` – Visualize a single point cloud in 3D with save-to-file support.
  - `plot_comparison()` – Side-by-side 3D comparison of Poisson vs Uniform sampling.
  - `plot_density_comparison()` – 2D density projection comparisons across XY, XZ, and YZ planes.
- Comparison visualization images in `images/` directory.
- `CHANGELOG.md` for tracking project history.

### Changed
- Updated README with sampling comparison visualizations, improved structure, and better documentation.
- Numba-optimized Poisson Disk Sampling (~140x speedup over pure Python).

## [0.1.0] - Initial Release

### Added
- Poisson Disk Sampling using Bridson's algorithm.
- Uniform random sampling in 3D bounding box.
- Basic 3D visualization with Matplotlib.
- Numba JIT compilation for performance.
