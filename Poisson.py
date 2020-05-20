#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 02:41:56 2020

@author: sumesh
@website : sumeshthakr.github.io
"""

import numpy as np
import cupy as cp
import pptk
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Poisson:
    
    def __init__(self, number_of_points, radius, max_x, max_y, max_z):
        
        self.k = number_of_points
        self.r = radius
        
        self.max_length = max_x
        self.max_width = max_y
        self.max_depth = max_z
        
        self.a = radius/np.sqrt(3)
        
        #Voxelization 
        #number of voxels
        
        self.nx, self.ny, self.nz = int(max_x/self.a)+1, int(max_y/self.a)+1, int(max_z/self.a)+1
        
        self.reset_voxel()
     
        
     #reseting vozels and initializing as null
     #Step 0
    def reset_voxel(self):
        
         
        cordinates_xyz = [(x_cor, y_cor, z_cor) for x_cor in range(self.nx)
                                                for y_cor in range(self.ny)
                                                for z_cor in range(self.nz)]
         
         #initiaizing to Null
         #Creating a dictionary of all the voxels
        self.voxels = {coords : None for coords in cordinates_xyz}
         
         
    def get_cell_coords(self, point):
        """Get the coordinates of the voxel in which point = (x,y,z) lies in."""

        return int(point[0] // self.a), int(point[1] // self.a), int(point[2] // self.a)
    
    def get_neighbour(self, coords):
        
        
        dxdydz = [(-1,-2,-1),(0,-2,-1),(1,-2, -1),(-2,-1, -1),(-1,-1, -1),(0,-1, -1),(1,-1, -1),(2,-1, -1),
                (-2,0, -1),(-1,0, -1),(1,0, -1),(2,0, -1),(-2,1,-1),(-1,1,-1),(0,1,-1),(1,1,-1),(2,1,-1),
                (-1,2,-1),(0,2,-1),(1,2,-1),(0,0,-1),(-1,-2, 0),(0,-2, 0),(1,-2, 0),(-2,-1, 0),(-1,-1, 0),(0,-1, 0),(1,-1, 0),(2,-1,0),
                (-2,0, 0),(-1,0,0),(1,0,0),(2,0,0),(-2,1,0),(-1,1,0),(0,1,0),(1,1,0),(2,1,0),
                (-1,2,0),(0,2,0),(1,2,0),(0,0,0),(-1,-2,1),(0,-2, 1),(1,-2, 1),(-2,-1, 1),(-1,-1,1),(0,-1,1),(1,-1,1),(2,-1,1),
                (-2,0,1),(-1,0,1),(1,0,1),(2,0,1),(-2,1,1),(-1,1,1),(0,1,1),(1,1,1),(2,1,1),
                (-1,2,1),(0,2,1),(1,2,1),(0,0,1)]
        
        neighbours = []
        for dx, dy, dz in dxdydz:
            neighbour_coords = coords[0] + dx, coords[1] + dy, coords[2] + dz
            if not (0 <= neighbour_coords[0] < self.nx and
                    0 <= neighbour_coords[1] < self.ny and
                    0 <= neighbour_coords[2] < self.nz) :
                # We're off the grid: no neighbours here.
                continue
            neighbour_cell = self.voxels[neighbour_coords]
            if neighbour_cell is not None:
                # This cell is occupied: store the index of the contained point
                neighbours.append(neighbour_cell)
        return neighbours
    
        
    def point_valid(self, pt):
        

        cell_coords = self.get_cell_coords(pt)
        for idx in self.get_neighbour(cell_coords):
            nearby_pt = self.samples[idx]
            # Squared distance between candidate point, pt, and this nearby_pt.
            distance2 = (nearby_pt[0]-pt[0])**2 + (nearby_pt[1]-pt[1])**2 + (nearby_pt[2]-pt[2])**2
            if distance2 < self.r**2:
                # The points are too close, so pt is not a candidate.
                return False
        # All points tested: if we're here, pt is valid
        return True
    
    def get_point(self, refpt):

        i = 0
        while i < self.k:
            rho, theta, yaw = (np.random.uniform(self.r, 2*self.r),
                          np.random.uniform(0, 2*np.pi),
                          np.random.uniform(0, 2*np.pi))
            pt = refpt[0] + rho*np.sin(theta)*np.cos(yaw), refpt[1] + rho*np.sin(theta)*np.sin(yaw) , refpt[2] + rho*np.cos(theta)
            if not (0 <= pt[0] < self.max_length and 0 <= pt[1] < self.max_width and 0<= pt[2] < self.max_depth):
                # This point falls outside the domain, so try again.
                continue
            if self.point_valid(pt):
                return pt
            i += 1
        # We failed to find a suitable point in the vicinity of refpt.
        return False
       
    def sample(self):
   
        pt = (np.random.uniform(0, self.max_length),
              np.random.uniform(0, self.max_width),
              np.random.uniform(0, self.max_depth))
        self.samples = [pt]
     
        self.voxels[self.get_cell_coords(pt)] = 0
   
        active = [0]

        while active:
      
            idx = np.random.choice(active)
            refpt = self.samples[idx]
            # Try to pick a new point relative to the reference point.
            pt = self.get_point(refpt)
            if pt:
                # Point pt is valid: add it to samples list and mark as active
                self.samples.append(pt)
                nsamples = len(self.samples) - 1
                active.append(nsamples)
                self.voxels[self.get_cell_coords(pt)] = nsamples
            else:
                # We had to give up looking for valid points near refpt, so
                # remove it from the list of "active" points.
                active.remove(idx)

        return self.samples

#Plotting

"""         
poisson = Poisson(30, 1.4, 10, 10, 10)    

sample = poisson.sample()

sample_np = np.array(sample)
x = sample_np[:,0]
y = sample_np[:,1]
z = sample_np[:,2]    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
ax.scatter(x, y, -z, zdir='z', c= 'red')
plt.show()"""
