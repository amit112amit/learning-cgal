#include <Eigen/Dense>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>

typedef Eigen::Vector3d Vector3d;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::Matrix3d Matrix3d;
typedef Eigen::Matrix3Xd Matrix3Xd;
typedef Eigen::Map<Matrix3Xd> Map3Xd;

int main(){
    clock_t t1;
    t1 = clock();
    for(auto i=0; i < 10000; ++i){
        auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName("T7.vtk");
        reader->Update();
        auto poly = reader->GetOutput();
        auto N = poly->GetNumberOfPoints();
        auto pts = (double*) poly->GetPoints()->GetData()->GetVoidPointer(0);
        Map3Xd points(pts,3,N);

        // Project points to unit sphere
        points.colwise().normalize();

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
        for( auto i=0; i < N-1; ++i ){
            proj.col(i) = ((p0(2) - l0(2,i))/l(2,i))*l.col(i) + l0.col(i);
            proj(i,2) = 0.0;
        }
        // Calculate the 2d delaunay triangulations
        auto pts2dArr = vtkSmartPointer<vtkDoubleArray>::New();
        pts2dArr->SetVoidArray((void*)proj.data(), 3*(N-1), 1);
        pts2dArr->SetNumberOfComponents(3);
        auto pts2d = vtkSmartPointer<vtkPoints>::New();
        pts2d->SetData(pts2dArr);
        auto poly2d = vtkSmartPointer<vtkPolyData>::New();
        poly2d->SetPoints(pts2d);
        auto d2d = vtkSmartPointer<vtkDelaunay2D>::New();
        d2d->SetInputData(poly2d);
        d2d->Update();
    }
    float diff((float)clock() - (float)t1);
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    //std::cout << dt.number_of_vertices() << std::endl;

    //std::ofstream out("mesh.off");
    //CGAL::export_triangulation_2_to_off( out, dt );
    return 0;
}
