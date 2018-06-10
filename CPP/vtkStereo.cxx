#include <Eigen/Dense>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkIdFilter.h>
#include <vtkPointData.h>

typedef Eigen::Vector3d Vector3d;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::Matrix3d Matrix3d;
typedef Eigen::Matrix3Xd Matrix3Xd;
typedef Eigen::Map<Matrix3Xd> Map3Xd;

int main(){
    clock_t t1;
    t1 = clock();
    for(auto i=0; i < 10000; ++i){
        vtkNew<vtkPolyDataReader> reader;
        reader->SetFileName("T7.vtk");
        reader->Update();
        auto poly = reader->GetOutput();
        auto N = poly->GetNumberOfPoints();
        auto pts = (double*) poly->GetPoints()->GetData()->GetVoidPointer(0);
        Map3Xd points(pts,3,N);

        // Project points to unit sphere
        points.colwise().normalize();

        // Reset the center of the sphere to origin by translating
        Vector3d center = points.rowwise().mean();
        points = points.colwise() - center;

        // Rotate all points so that the point in 0th column is along z-axis
        Vector3d c = points.col(0);
        double_t cos_t = c(2);
        double_t sin_t = std::sqrt( 1 - cos_t*cos_t );
        Vector3d axis;
        axis << c(1), -c(0), 0.;
        Matrix3d rotMat, axis_cross, outer;
        axis_cross << 0. , -axis(2), axis(1),
                        axis(2), 0., -axis(0),
                        -axis(1), axis(0), 0.;

        outer.noalias() = axis*axis.transpose();

        rotMat = cos_t*Matrix3d::Identity() + sin_t*axis_cross + (1-cos_t)*outer;
        Matrix3Xd rPts(3,N);
        rPts = rotMat*points; // The points on a sphere rotated

        // Calculate the stereographic projections
        Vector3d p0;
        Map3Xd l0( &(rPts(0,1)), 3, N-1 );
        Matrix3Xd l(3,N-1), proj(3,N-1);
        p0 << 0,0,-1;
        c = rPts.col(0);
        l = (l0.colwise() - c).colwise().normalized();
        for( auto j=0; j < N-1; ++j ){
            proj.col(j) = ((p0(2) - l0(2,j))/l(2,j))*l.col(j) + l0.col(j);
            proj(j,2) = 0.0;
        }
        // Calculate the 2d delaunay triangulations
        vtkNew<vtkDoubleArray> pts2dArr;
        pts2dArr->SetVoidArray((void*)proj.data(), 3*(N-1), 1);
        pts2dArr->SetNumberOfComponents(3);
        vtkNew<vtkPoints> pts2d;
        pts2d->SetData(pts2dArr);
        vtkNew<vtkPolyData> poly2d;
        poly2d->SetPoints(pts2d);
        vtkNew<vtkIdFilter> idf;
        idf->PointIdsOn();
        idf->SetIdsArrayName("OrigIds");
        idf->SetInputData(poly2d);
        vtkNew<vtkDelaunay2D> d2d;
        d2d->SetInputConnection(idf->GetOutputPort());
        d2d->Update();

        // Write the triangulation to file
        vtkNew<vtkCellArray> final;
        auto stereoTris = d2d->GetOutput()->GetPolys();
        auto idArr = d2d->GetOutput()->GetPointData()->GetArray("OrigIds");
        vtkNew<vtkIdList> idL;
        stereoTris->InitTraversal();
        while( stereoTris->GetNextCell(idL) ){
            final->InsertNextCell(3);
            for(auto j=0; j < idL->GetNumberOfIds(); ++j)
                final->InsertCellPoint( int(
                                        idArr->GetTuple1(idL->GetId(j))) );
        }

    }
    float diff((float)clock() - (float)t1);
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    return 0;
}
