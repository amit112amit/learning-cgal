#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <Eigen/Dense>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Surface_mesh<Point> Surface_mesh;
typedef Surface_mesh::Vertex_index Vertex;
typedef Surface_mesh::Edge_index Edge;
typedef Surface_mesh::Face_index Face;
typedef Surface_mesh::Halfedge_index Halfedge;
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
        reader->SetFileName("Bad.vtk");
        reader->Update();
        auto poly = reader->GetOutput();
        auto N = poly->GetNumberOfPoints();
        auto pts = static_cast<double_t*>(poly->GetPoints()->GetData()->GetVoidPointer(0));
        Map3Xd points(pts,3,N);

        // Project points to unit sphere
        points.colwise().normalize();

        // Reset the center of the sphere to origin by translating
        Vector3d center = points.rowwise().mean();
        points = points.colwise() - center;

        std::vector<Point> spherePoints(N);
        for( auto i=0; i < N; ++i){
            spherePoints[i] = Point(points(0,i),points(1,i),points(2,i));
        }

        // Make a mesh object
        Surface_mesh sm;

        // Calculate the convex hull
        CGAL::convex_hull_3(spherePoints.begin(),
                            spherePoints.end(), sm);

        // To extract the surface
        // Write to a vtk file
        vtkNew<vtkCellArray> triangles;
        for( auto f : sm.faces() ){
            triangles->InsertNextCell(3);
	    auto he = sm.halfedge( f );
            for( auto j=0; j < 3; ++j){
                triangles->InsertCellPoint( static_cast<size_t>(sm.target(he)) );
		he = sm.next( he );
            }
        }
        poly->SetPolys(triangles);
        vtkNew<vtkPolyDataWriter> writer;
        writer->SetFileName("CH_Mesh.vtk");
        writer->SetInputData(poly);
        writer->Write();
    }
    float diff(static_cast<float>(clock()) - static_cast<float>(t1));
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    return 0;
}
