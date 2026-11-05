#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 10:12:44 2018

@author: cecidip
"""
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def load_vel_vtr(vtr_path):
    """
    Reads the LaRonde .vtr file and extracts Vp as a flat 1D array of nodes.
    Matches the exact output structure of the SGS script.
    """
    reader = vtk.vtkXMLRectilinearGridReader()
    reader.SetFileName(vtr_path)
    reader.Update()
    grid = reader.GetOutput()
    
    # Extract Vp from PointData (nodes) and cast to float32
    vp_vtk = grid.GetPointData().GetArray("Vp")
    vp_flat = vtk_to_numpy(vp_vtk).astype(np.float32)
    
    return vp_flat

def plot_rnSGS(g, E):
    """
    Plots slices of randomly chosen ensemble members.
    If you want to save the plot call the function as:
    fig = plot_rnSGS(g, E)
    fig.savefig("somefile.png")
    """
    x = g.x
    y = g.y
    z = g.z

    N = len(E) 
    # Grab random ensemble members and reshape to 3D grid
    Vrn1 = E[np.random.randint(N), :].reshape(g.shape)
    Vrn2 = E[np.random.randint(N), :].reshape(g.shape)
    Vrn3 = E[np.random.randint(N), :].reshape(g.shape)
    Vrn4 = E[np.random.randint(N), :].reshape(g.shape)
    
    fig = plt.figure(figsize=(10, 8))
    
    # Swapped pcolor for pcolormesh (much faster and standard in modern Matplotlib)
    plt.subplot(221)
    plt.pcolormesh(x, z, np.squeeze(Vrn1[:, 10, :].T), cmap='jet', shading='auto')
    plt.gca().invert_yaxis()
    plt.xlabel('X')
    plt.ylabel('Z')
    plt.colorbar()
    
    plt.subplot(222)
    plt.pcolormesh(x, z, np.squeeze(Vrn2[:, 10, :].T), cmap='jet', shading='auto')
    plt.gca().invert_yaxis()
    plt.xlabel('Y')
    plt.ylabel('Z')
    plt.colorbar()
    
    plt.subplot(223)
    plt.pcolormesh(x, z, np.squeeze(Vrn3[:, 10, :].T), cmap='jet', shading='auto')
    plt.gca().invert_yaxis()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar()
    
    plt.subplot(224)
    plt.pcolormesh(x, z, np.squeeze(Vrn4[:, 10, :].T), cmap='jet', shading='auto')
    plt.gca().invert_yaxis()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar()
    
    plt.suptitle('Random SGS Realizations')
    plt.tight_layout()
    plt.show()
    return fig

'''
def esfera(g, centro, V, pp=0.05, rr=120, vlim=None):
    """
    Applies a spherical velocity perturbation around the hypocenter.
    
    """
    x0, y0, z0 = centro[0], centro[1], centro[2]
    
    # Calculate the center points of the cells
    dx = g.x[1] - g.x[0]
    dy = g.y[1] - g.y[0]
    dz = g.z[1] - g.z[0]
    
    xc = g.x[:-1] + dx / 2.0
    yc = g.y[:-1] + dy / 2.0
    zc = g.z[:-1] + dz / 2.0
    
    # Generate 3D coordinate matrices based on cell centers
    X, Y, Z = np.meshgrid(xc, yc, zc, indexing='ij')
    
    # Calculate distance for the entire grid at once
    dis = np.sqrt((X - x0)**2 + (Y - y0)**2 + (Z - z0)**2)
    
    # Apply perturbation percentage where distance is within the radius
    corr = np.ones(g.shape)
    corr[dis <= rr] = (1 + pp)
    
    corrV = corr.flatten()
    V_updated = V * corrV
    
    if vlim is not None:
        V_updated = np.clip(V_updated, a_min=vlim, a_max=None)
        
    return V_updated
'''
@njit(cache=True)
def esfera(xc, yc, zc, V, x0, y0, z0, pp, rr):
    nx, ny, nz = len(xc), len(yc), len(zc)
    V_updated = np.empty(nx * ny * nz, dtype=V.dtype)
    idx = 0
    rr_sq = rr ** 2 # Avoids expensive square roots during the loop
    
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                dist_sq = (xc[i] - x0)**2 + (yc[j] - y0)**2 + (zc[k] - z0)**2
                if dist_sq <= rr_sq:
                    V_updated[idx] = V[idx] * (1.0 + pp)
                else:
                    V_updated[idx] = V[idx]
                idx += 1
                
    return V_updated


def stats(E, Vtrue):
    """
    Computes ensemble statistics against the true velocity model.
    Reference: Tarrahi 2015 "Integration of microseismic monitoring data.. Kalman filter"
    """
    N, m = E.shape 
    # RMSE calculation
    tmp = np.sum((E - Vtrue)**2, axis=0)
    RMSE = np.sum(np.sqrt(tmp / N)) / m 
    # Standard deviation calculation
    vmean = np.sum(E, axis=0) / N
    tmp1 = np.sum((E - vmean)**2, axis=0)
    Sp = np.sum(np.sqrt(tmp1 / N)) / m
    
    return RMSE, Sp