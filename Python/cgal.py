import warnings
warnings.filterwarnings('ignore')
import numpy as np
import vtk as v
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

# Find the angle of sphere[0] with z-axis and rotate all points
c = sphere[0]
cos_t = c[2]
sin_t = np.sqrt(1 - cos_t**2)
axis = [c[1],-c[0],0.]

# The rotation matrix
axis_cross = np.array([[0., -axis[2], axis[1]],
              [axis[2], 0., -axis[0]],
              [-axis[1], axis[0], 0.]])
rotMat = cos_t*np.eye(3) + sin_t*axis_cross + (1 - cos_t)*np.outer(axis,axis)

# Rotate the points and make a vtk array
rPts = np.matmul(rotMat,sphere.T).T
rPd = numpy2PolyData( rPts )
writeVTK( rPd, 'Rotated.vtk')

# Now we will calculate the stereographic projection
p0 = np.array([[0.,0.,-2]])
c = rPts[0]
l0 = rPts[1:,:]
l = (l0 - c)/np.linalg.norm( l0-c, axis=1)[:,np.newaxis]
d = (p0[:,2] - l0[:,2])/l[:,2]
proj3 = d[:,np.newaxis]*l + l0
proj = proj3.copy()

# Write the stereographic projection set up to vtk
proj[:,2] = 0
sPd = numpy2PolyData( proj )
idf = v.vtkIdFilter()
idf.PointIdsOn()
idf.SetIdsArrayName('OrigIds')
idf.SetInputData( sPd )
d2d = v.vtkDelaunay2D()
d2d.SetInputConnection( idf.GetOutputPort() )
d2d.Update()
meshTemp = d2d.GetOutput().GetPolys()
meshTemp.InitTraversal()
usg = v.vtkUnstructuredGrid()
usg.Allocate(198,198)
ids = v.vtkIdList()
origIds = d2d.GetOutput().GetPointData().GetArray('OrigIds')
while( meshTemp.GetNextCell( ids ) ):
    n = ids.GetNumberOfIds()
    uIds = v.vtkIdList()
    for i in range(n):
        uIds.InsertNextId( int(origIds.GetTuple1(ids.GetId(i))) )
    usg.InsertNextCell(5,uIds)

# Insert Projection Lines
uPts = v.vtkPoints()
uPtsArr = np.vstack( [proj3,rPts] )
uPtsArr[N-1] = c
uPts.SetData( dsa.numpyTovtkDataArray( uPtsArr ) )
for i in range(N-1):
    uIds = v.vtkIdList()
    uIds.InsertNextId( N-1 )
    uIds.InsertNextId( i )
    uIds.InsertNextId( i+N )
    usg.InsertNextCell( 4, uIds )

usg.SetPoints( uPts )
wr = v.vtkUnstructuredGridWriter()
wr.SetFileName('Scene.vtk')
wr.SetInputData(usg)
wr.Write()

# Now write the final spherical mesh to polydate file

