#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 19:20:28 2020

@author: sumesh
"""

"""Compatibility wrapper for uniform sampling script."""
from pointcloud_sampling import sample_uniform

if __name__ == "__main__":
    pts = sample_uniform(100, (0, 10), (0, 10), (0, 10))
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
