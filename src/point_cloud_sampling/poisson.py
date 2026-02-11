import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import List, Tuple, Optional, Dict
from numba import njit
import math

@njit
def get_cell_coords(point, a):
    return (int(point[0] // a), int(point[1] // a), int(point[2] // a))

@njit
def generate_poisson_points(max_attempts, radius, max_dims):
    # Unpack max dimensions
    max_x, max_y, max_z = max_dims
    r2 = radius * radius
    a = radius / math.sqrt(3)
    
    # Grid dimensions
    grid_w = int(math.ceil(max_x / a))
    grid_h = int(math.ceil(max_y / a))
    grid_d = int(math.ceil(max_z / a))
    
    # Dense grid: -1 indicates empty
    # flattened index = x + y * w + z * w * h
    grid_size = grid_w * grid_h * grid_d
    grid = np.full(grid_size, -1, dtype=np.int32)
    
    # Preallocate arrays
    # Estimate max points: Volume / (minimal spacing volume per point approx)
    # A generous upper bound is Volume / (r/2)^3 or similar.
    # Volume = max_x * max_y * max_z
    # Each point occupies roughly (r)^3 space loosely.
    # Let's be safe and allocate plenty, or use dynamic resizing.
    # Dynamic resizing is possible but simpler to just allocate a large buffer.
    # 50x50x50 / 1.5^3 ~ 37000. Let's alloc 10x that to be safe or just 1M.
    # Since we return a sliced array, memory usage is temporary.
    max_points = int((max_x * max_y * max_z) / (a * a * a) * 2) + 1000
    if max_points < 1000: max_points = 1000
    
    samples = np.empty((max_points, 3), dtype=np.float64)
    sample_count = 0
    
    active_indices = np.empty(max_points, dtype=np.int32)
    active_count = 0
    
    # Step 1: Initialize with a random point
    # Use array for point to avoid object creation
    first_pt_x = np.random.uniform(0, max_x)
    first_pt_y = np.random.uniform(0, max_y)
    first_pt_z = np.random.uniform(0, max_z)
    
    samples[0, 0] = first_pt_x
    samples[0, 1] = first_pt_y
    samples[0, 2] = first_pt_z
    sample_count += 1
    
    active_indices[0] = 0
    active_count += 1
    
    cx = int(first_pt_x // a)
    cy = int(first_pt_y // a)
    cz = int(first_pt_z // a)
    
    grid_idx = cx + cy * grid_w + cz * grid_w * grid_h
    grid[grid_idx] = 0
    
    while active_count > 0:
        rand_idx = np.random.randint(0, active_count)
        pt_idx = active_indices[rand_idx]
        
        ref_px = samples[pt_idx, 0]
        ref_py = samples[pt_idx, 1]
        ref_pz = samples[pt_idx, 2]
        
        found = False
        for _ in range(max_attempts):
            # Generate random point around ref_pt
            rho = np.random.uniform(radius, 2 * radius)
            theta = np.random.uniform(0, 2 * math.pi)
            phi = np.random.uniform(0, math.pi)
            
            x = rho * math.sin(phi) * math.cos(theta)
            y = rho * math.sin(phi) * math.sin(theta)
            z = rho * math.cos(phi)
            
            cand_x = ref_px + x
            cand_y = ref_py + y
            cand_z = ref_pz + z
            
            # Check bounds
            if not (0 <= cand_x < max_x and 
                    0 <= cand_y < max_y and 
                    0 <= cand_z < max_z):
                continue
            
            # Check neighbors
            cx = int(cand_x // a)
            cy = int(cand_y // a)
            cz = int(cand_z // a)
            
            # Check 5x5x5 neighborhood around cell
            valid = True
            
            # Optimized range loop
            z_min = max(0, cz - 2)
            z_max = min(grid_d, cz + 3)
            y_min = max(0, cy - 2)
            y_max = min(grid_h, cy + 3)
            x_min = max(0, cx - 2)
            x_max = min(grid_w, cx + 3)
            
            for nz in range(z_min, z_max):
                for ny in range(y_min, y_max):
                    for nx in range(x_min, x_max):
                        idx = nx + ny * grid_w + nz * grid_w * grid_h
                        neighbor_pt_idx = grid[idx]
                        
                        if neighbor_pt_idx != -1:
                            # Distance check
                            # Squared distance
                            dx = cand_x - samples[neighbor_pt_idx, 0]
                            dy = cand_y - samples[neighbor_pt_idx, 1]
                            dz = cand_z - samples[neighbor_pt_idx, 2]
                            dist_sq = dx*dx + dy*dy + dz*dz
                            
                            if dist_sq < r2:
                                valid = False
                                break # break x loop
                    if not valid: break # break y loop
                if not valid: break # break z loop
            
            if valid:
                if sample_count >= max_points:
                    # buffer full, this is an error condition for this simple version
                    # In real production code, we'd resize.
                    # For now just break or return what we have? 
                    # Or raise error (but njit doesn't like exceptions much)
                    # Let's just stop adding points to be safe.
                    found = False # treat as not found to exit? No, just break main loop
                    active_count = 0 # force exit
                    break

                samples[sample_count, 0] = cand_x
                samples[sample_count, 1] = cand_y
                samples[sample_count, 2] = cand_z
                
                active_indices[active_count] = sample_count
                
                grid_idx = cx + cy * grid_w + cz * grid_w * grid_h
                grid[grid_idx] = sample_count
                
                sample_count += 1
                active_count += 1
                
                found = True
                break
        
        if not found and active_count > 0:
            # Remove from active list (swap and pop)
            active_indices[rand_idx] = active_indices[active_count - 1]
            active_count -= 1
            
    return samples[:sample_count]

class PoissonSampler:
    """
    Poisson Disk Sampling in 3D (Optimized with Numba).
    """
    
    def __init__(self, max_attempts: int, radius: float, max_x: float, max_y: float, max_z: float):
        self.k = max_attempts
        self.r = radius
        self.max_dims = np.array([max_x, max_y, max_z], dtype=np.float64)

    def sample(self) -> np.ndarray:
        """Generate samples."""
        return generate_poisson_points(self.k, self.r, self.max_dims)

def plot_samples(samples: np.ndarray):
    """Visualize samples."""
    if samples.size == 0:
        print("No samples to plot.")
        return
        
    x = samples[:, 0]
    y = samples[:, 1]
    z = samples[:, 2]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.scatter(x, y, z, c='r', marker='.')
    plt.show()

if __name__ == "__main__":
    # Example usage
    sampler = PoissonSampler(30, 1.4, 10, 10, 10)
    points = sampler.sample()
    print(f"Generated {len(points)} points.")
    # plot_samples(points)
