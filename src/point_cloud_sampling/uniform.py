import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import Optional, Tuple

class UniformSampler:
    """
    Uniform Sampler in 3D.
    """
    
    @staticmethod
    def sample(n_points: int, x_range: Tuple[float, float], y_range: Tuple[float, float], z_range: Tuple[float, float]) -> np.ndarray:
        """
        Generate n_points uniformly distributed in the defined box.
        
        Args:
            n_points (int): Number of points to generate.
            x_range (Tuple[float, float]): (min, max) for X.
            y_range (Tuple[float, float]): (min, max) for Y.
            z_range (Tuple[float, float]): (min, max) for Z.
            
        Returns:
            np.ndarray: Array of shape (n_points, 3) containing coordinates.
        """
        xs = np.random.uniform(x_range[0], x_range[1], n_points)
        ys = np.random.uniform(y_range[0], y_range[1], n_points)
        zs = np.random.uniform(z_range[0], z_range[1], n_points)
        
        return np.column_stack((xs, ys, zs))

def plot_samples(samples: np.ndarray):
    """Visualize samples."""
    if samples.size == 0:
        print("No samples to plot.")
        return

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.scatter(samples[:, 0], samples[:, 1], samples[:, 2], marker='.')
    plt.show()

if __name__ == "__main__":
    # Example usage
    n = 100
    points = UniformSampler.sample(n, (0, 10), (0, 10), (-50, -25))
    print(f"Generated {len(points)} points.")
    plot_samples(points)