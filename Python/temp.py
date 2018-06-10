import vtk as v
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.numpy_support import numpy_to_vtk
from vtk.numpy_interface import algorithms as alg
from scipy.spatial import Delaunay

# Get input and output polydata
reader = v.vtkPolyDataReader()
reader.SetFileName('T7.vtk')
reader.Update()
pdi = reader.GetOutput()
pdo = v.vtkPolyData()

# Copy points from input to output
pdo.SetPoints( pdi.GetPoints() )

# Project the points to a sphere
pd = dsa.WrapDataObject(pdi)
sphere = alg.norm( pd.Points )

# Reset the center of sphere to origin
sphere -= np.mean(sphere,axis=0)

# Calculate the spherical angles of the points
theta = np.arctan2(sphere[:,1],sphere[:,0])
phi = np.arctan2(np.sqrt(sphere[:,0]**2 + sphere[:,1]**2),sphere[:,2])

# Delaunay triangulation using SciPy
conn = Delaunay(np.hstack([theta[:,np.newaxis],phi[:,np.newaxis]]))

# Set the connectivity matrix to output polydata
triangles = v.vtkCellArray()
triData = np.hstack( [3*np.ones((conn.simplices.shape[0],1)),conn.simplices] )
vtkArr = numpy_to_vtk( triData, deep=1, array_type=v.VTK_ID_TYPE)
triangles.SetCells( triData.shape[0], vtkArr )
pdo.SetPolys( triangles )
