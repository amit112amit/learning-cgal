#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <Eigen/Dense>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned,K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Delaunay;
typedef Delaunay::Face_circulator Face_circulator;
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
        reader->SetFileName("Bad.vtk");
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
        }
        // Insert the projected points in a CGAL vertex_with_info vector
        std::vector< std::pair< Point, unsigned> > verts;
        for( auto j=0; j < N-1; ++j ){
            verts.push_back(std::make_pair(Point(proj(0,j),proj(1,j)),j+1));
        }

        Delaunay dt( verts.begin(), verts.end() );

        // Write the finite faces of the triangulation to a VTK file
        vtkNew<vtkCellArray> triangles;
        for( auto ffi = dt.finite_faces_begin(); ffi != dt.finite_faces_end(); ++ffi){
            triangles->InsertNextCell(3);
            for(auto j=2; j >= 0; --j)
                triangles->InsertCellPoint(ffi->vertex(j)->info());
        }

        // Iterate over infinite faces
        Face_circulator fc = dt.incident_faces(dt.infinite_vertex()), done(fc);
        if (fc != 0) {
            do{
                triangles->InsertNextCell(3);
                for(auto j=2; j >= 0; --j){
                    auto vh = fc->vertex(j);
                    auto id = dt.is_infinite(vh)? 0 : vh->info();
                    triangles->InsertCellPoint(id);
                }
            }while(++fc != done);
        }
        poly->SetPolys(triangles);

        // Write to VTK file
        vtkNew<vtkPolyDataWriter> writer;
        writer->SetFileName("StereoMesh.vtk");
        writer->SetInputData(poly);
        writer->Write();
    }
    float diff((float)clock() - (float)t1);
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    return 0;
}
