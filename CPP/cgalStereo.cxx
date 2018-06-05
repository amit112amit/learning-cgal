#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Triangulation_off_ostream_2.h>
#include <Eigen/Dense>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3 Point;
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

        /*
    // Write the rotated points to file
    auto rPtsArr = vtkSmartPointer<vtkDoubleArray>::New();
    auto rPtsPts = vtkSmartPointer<vtkPoints>::New();
    auto rPtsPd = vtkSmartPointer<vtkPolyData>::New();
    auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();

    rPtsArr->SetVoidArray( (void*)rPts.data(), 3*N, 1 );
    rPtsArr->SetNumberOfComponents(3);
    rPtsPts->SetData( rPtsArr );
    rPtsPd->SetPoints( rPtsPts );
    writer->SetFileName("Rotated.vtk");
    writer->SetInputData( rPtsPd );
    writer->Write();
    */

        // Calculate the stereographic projections
        Vector3d p0;
        Map3Xd l0( &(rPts(0,1)), 3, N-1 );
        Matrix3Xd l(3,N-1), proj(3,N-1);
        p0 << 0,0,-1;
        c = rPts.col(0);
        l = (l0.colwise() - c).colwise().normalized();
        for( auto i=0; i < N-1; ++i ){
            proj.col(i) = ((p0(2) - l0(2,i))/l(2,i))*l.col(i) + l0.col(i);
        }
        // Insert the projected points in a CGAL Point_3 vector
        std::vector< Point > verts(N-1);
        for( auto i=0; i < N-1; ++i ){
            verts[i] = Point( proj(0,i), proj(1,i), proj(2,i) );
        }

        Delaunay dt( verts.begin(), verts.end() );
    }
    float diff((float)clock() - (float)t1);
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    //std::cout << dt.number_of_vertices() << std::endl;

    //std::ofstream out("mesh.off");
    //CGAL::export_triangulation_2_to_off( out, dt );
    return 0;
}
