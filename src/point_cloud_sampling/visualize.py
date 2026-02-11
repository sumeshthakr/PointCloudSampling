import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import Optional


def plot_samples(samples: np.ndarray, title: str = "Point Cloud",
                 color: str = 'steelblue', marker_size: float = 2.0,
                 save_path: Optional[str] = None):
    """Visualize a single point cloud in 3D.

    Args:
        samples: Array of shape (N, 3) with point coordinates.
        title: Plot title.
        color: Marker color.
        marker_size: Marker size.
        save_path: If provided, save the figure to this path instead of showing.
    """
    if samples.size == 0:
        print("No samples to plot.")
        return

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(samples[:, 0], samples[:, 1], samples[:, 2],
               c=color, marker='.', s=marker_size, alpha=0.6)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()


def plot_comparison(poisson_samples: np.ndarray, uniform_samples: np.ndarray,
                    save_path: Optional[str] = None):
    """Create a side-by-side comparison of Poisson and Uniform sampling.

    Args:
        poisson_samples: Array of shape (N, 3) from Poisson sampling.
        uniform_samples: Array of shape (M, 3) from Uniform sampling.
        save_path: If provided, save the figure to this path instead of showing.
    """
    fig = plt.figure(figsize=(16, 6))

    # Poisson Disk Sampling
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.scatter(poisson_samples[:, 0], poisson_samples[:, 1],
                poisson_samples[:, 2], c='steelblue', marker='.', s=2, alpha=0.6)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.set_title(f'Poisson Disk Sampling\n({len(poisson_samples)} points)')

    # Uniform Sampling
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.scatter(uniform_samples[:, 0], uniform_samples[:, 1],
                uniform_samples[:, 2], c='coral', marker='.', s=2, alpha=0.6)
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Z')
    ax2.set_title(f'Uniform Sampling\n({len(uniform_samples)} points)')

    plt.suptitle('Sampling Method Comparison', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()


def plot_density_comparison(poisson_samples: np.ndarray, uniform_samples: np.ndarray,
                            save_path: Optional[str] = None):
    """Create 2D density projection comparisons (XY, XZ, YZ planes).

    Projects 3D points onto 2D planes to visualize point distribution density.

    Args:
        poisson_samples: Array of shape (N, 3) from Poisson sampling.
        uniform_samples: Array of shape (M, 3) from Uniform sampling.
        save_path: If provided, save the figure to this path instead of showing.
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    planes = [('X', 'Y', 0, 1), ('X', 'Z', 0, 2), ('Y', 'Z', 1, 2)]

    for col, (xlabel, ylabel, xi, yi) in enumerate(planes):
        # Poisson
        axes[0, col].scatter(poisson_samples[:, xi], poisson_samples[:, yi],
                             c='steelblue', marker='.', s=1, alpha=0.4)
        axes[0, col].set_xlabel(xlabel)
        axes[0, col].set_ylabel(ylabel)
        axes[0, col].set_aspect('equal')
        if col == 0:
            axes[0, col].set_ylabel(f'Poisson\n{ylabel}', fontsize=12, fontweight='bold')
        axes[0, col].set_title(f'{xlabel}-{ylabel} Projection')

        # Uniform
        axes[1, col].scatter(uniform_samples[:, xi], uniform_samples[:, yi],
                             c='coral', marker='.', s=1, alpha=0.4)
        axes[1, col].set_xlabel(xlabel)
        axes[1, col].set_ylabel(ylabel)
        axes[1, col].set_aspect('equal')
        if col == 0:
            axes[1, col].set_ylabel(f'Uniform\n{ylabel}', fontsize=12, fontweight='bold')

    plt.suptitle('2D Density Projections: Poisson vs Uniform Sampling',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()
