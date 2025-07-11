import numpy as np
import scipy.io
import pyvista as pv

# Load MATLAB data
data = scipy.io.loadmat('results/liutex_data.mat')
x = data['x_cc'].squeeze()
y = data['y_cc'].squeeze()
z = data['z_cc'].squeeze()
vel1 = data['vel1']
pres = data['pres']
liutex = data['liutex']

mp, np_, pp = vel1.shape

# Create meshgrid of coordinates (note: MATLAB order is (x,y,z) = (i,j,k))
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')  # shape: (mp, np, pp)

# Reshape to 1D arrays for VTK
points = np.column_stack([X.ravel(order='F'),
                          Y.ravel(order='F'),
                          Z.ravel(order='F')])

# Create structured grid
grid = pv.StructuredGrid()
grid.dimensions = (mp, np_, pp)
grid.points = points

# Add variables (flattened in Fortran order)
grid['vel1'] = vel1.ravel(order='F')
grid['pres'] = pres.ravel(order='F')
grid['liutex'] = liutex.ravel(order='F')

# Save to .vts file
grid.save('results/liutex_data.vts')
