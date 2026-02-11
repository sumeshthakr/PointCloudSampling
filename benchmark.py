import time
import numpy as np
from src.point_cloud_sampling.poisson import PoissonSampler

def benchmark(n_runs=5):
    print(f"Benchmarking Poisson Sampler over {n_runs} runs...")
    
    times = []
    total_points = 0
    
    # Parameters for a non-trivial workload
    # 50x50x50 box ~ 125,000 unit cubes. Radius 1.5.
    max_dims = (50, 50, 50)
    radius = 1.5
    
    for i in range(n_runs):
        start_time = time.time()
        sampler = PoissonSampler(30, radius, *max_dims)
        points = sampler.sample()
        end_time = time.time()
        
        duration = end_time - start_time
        times.append(duration)
        total_points += len(points)
        print(f"Run {i+1}: {duration:.4f}s ({len(points)} points)")
        
    avg_time = np.mean(times)
    avg_points = total_points / n_runs
    
    print(f"\nAverage Time: {avg_time:.4f}s")
    print(f"Average Points: {avg_points:.1f}")
    print(f"Points per Second: {avg_points / avg_time:.1f}")

if __name__ == "__main__":
    benchmark()
