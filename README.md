# PointCloudSampling

Python package providing Poisson-disk and uniform sampling for 3D point clouds.

## Requirements
- Python 3.8+ (virtualenv recommended)

## Install (editable / from source)
# From project root:
python -m venv .venv311
source .venv311/bin/activate
pip install -e .[vis,test]

## Usage
from pointcloud_sampling import Poisson, sample_uniform

# Poisson sampling
p = Poisson(k=30, r=1.4, max_x=10, max_y=10, max_z=10)
points = p.sample()  # (N,3) numpy array

# Uniform sampling in a box
pts = sample_uniform(100, (0,10), (0,10), (0,10))

## Tests
pip install -e .[test]
pytest

## Notes
- For visualization install the `vis` extra (adds `open3d`)
- This project is packaged via `pyproject.toml` for pip/setuptools.

---

## Publishing to PyPI (automated)
This repo includes a GitHub Actions workflow that will build and publish the package to PyPI when a tag matching `v*.*.*` is pushed.

1. Create a PyPI API token (on https://pypi.org/manage/account/token/).
2. Add the token as a GitHub repo secret named `PYPI_API_TOKEN` in the repository settings.
3. Create a tag and push it to GitHub, for example:

   git tag -a v0.1.0 -m "Release v0.1.0"
   git push origin v0.1.0

The workflow will run and publish the package automatically.

If you prefer, I can publish now for you — I'll need a PyPI API token (or you can set it in GitHub secrets and I'll create the tag).