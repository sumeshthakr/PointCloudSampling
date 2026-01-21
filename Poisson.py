"""Compatibility wrapper: import the packaged implementation and provide CLI."""
from pointcloud_sampling import Poisson


if __name__ == "__main__":
    poisson = Poisson(k=30, r=1.4, max_x=10, max_y=10, max_z=10)
    pts = poisson.sample()
    try:
        import open3d as o3d
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(pts)
        o3d.visualization.draw_geometries([pcd])
    except Exception:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2])
        plt.show()
