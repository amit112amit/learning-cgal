import warnings
warnings.filterwarnings('ignore')
import numpy as np
import matplotlib.pyplot as plt
import vtk as v
from scipy.spatial import Delaunay
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as alg

def numpy2PolyData(data):
    vtkArr = dsa.numpyTovtkDataArray( data )
    pts = v.vtkPoints()
    pts.SetData(vtkArr)
    pd = v.vtkPolyData()
    pd.SetPoints( pts )
    vgf = v.vtkVertexGlyphFilter()
    vgf.SetInputData( pd )
    vgf.Update()
    return vgf.GetOutput()

def writeVTK( pd, fileName ):
    wr = v.vtkDataSetWriter()
    wr.SetFileName(fileName)
    wr.SetInputData(pd)
    wr.Write()

def triangulate(sphereXyz):
    """
    Generates a triangle mesh for a spherical point cloud.
    """
    sphereXyz = np.around(sphereXyz,decimals=2)
    sphereArr = dsa.numpyTovtkDataArray( sphereXyz, name='SpherePts' )
    pts = v.vtkPoints()
    pts.SetData( sphereArr )
    sphere = v.vtkPolyData()
    sphere.SetPoints( pts )

    # Store the original point ids
    idf = v.vtkIdFilter()
    idf.SetIdsArrayName('PointIds')
    idf.PointIdsOn()
    idf.SetInputData( sphere )

    # Delaunay3D to make a convex hull
    d3d = v.vtkDelaunay3D()
    d3d.SetInputConnection( idf.GetOutputPort() )

    # Extract the surface
    surf = v.vtkDataSetSurfaceFilter()
    surf.SetInputConnection( d3d.GetOutputPort() )
    surf.Update()

    # Now make a new cell array mapping to the old ids
    polyCells = v.vtkCellArray()
    sphereCells = surf.GetOutput().GetPolys()
    sphereCells.InitTraversal()
    origIds = surf.GetOutput().GetPointData().GetArray('PointIds')
    ptIds = v.vtkIdList()
    while( sphereCells.GetNextCell( ptIds ) ):
        polyCells.InsertNextCell(3)
        polyCells.InsertCellPoint( int(origIds.GetTuple1( ptIds.GetId(0) )) )
        polyCells.InsertCellPoint( int(origIds.GetTuple1( ptIds.GetId(1) )) )
        polyCells.InsertCellPoint( int(origIds.GetTuple1( ptIds.GetId(2) )) )

    connectivity = dsa.vtkDataArrayToVTKArray( polyCells.GetData() )
    return connectivity

# Read T=7 points
rd = v.vtkPolyDataReader()
rd.SetFileName('T7.vtk')
rd.Update()
pdVtk = rd.GetOutput()
pd = dsa.WrapDataObject( pdVtk )
pts = pd.Points
N = pts.shape[0]

# Get the points on a sphere
sphere = alg.norm( pts )

# Set the center of the sphere to origin
center = np.mean(sphere,axis=0)
sphere -= center

conn = triangulate(sphere)

# Convert all points to spherical coordinates
theta = np.arctan2(sphere[:,1],sphere[:,0])
phi = np.arctan2(np.sqrt(sphere[:,0]**2 + sphere[:,1]**2),sphere[:,2])
print('Theta =')
print(theta)
print('Phi =')
print(phi)
phi = np.append(phi, phi)
theta = np.append(theta, theta - 2*np.pi)

conn2 = Delaunay(np.hstack([theta[:,np.newaxis],phi[:,np.newaxis]]))

fig, (ax1,ax2) = plt.subplots(2,1,sharex=True)
ax1.plot(theta,phi,'.')
ax1.triplot(theta,phi,triangles=conn2.simplices)
ax1.axvline(x=-np.pi)
ax1.set_ylabel(r'$\phi$')
ax2.plot(theta,phi,'.')
ax2.triplot(theta,phi,triangles=conn)
ax2.axvline(x=-np.pi)
ax2.set_xlabel(r'$\theta$')
ax2.set_ylabel(r'$\phi$')
plt.tight_layout()
plt.show()
