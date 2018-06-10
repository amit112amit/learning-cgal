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
        auto pts = (double*) poly->GetPoints()->GetData()->GetVoidPointer(0);
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

        /*
        // To extract the surface
        std::vector<Cell_handle> cells;
        T.tds().incident_cells( T.infinite_vertex(),
                                std::back_inserter(cells) );

        // Write to a vtk file
        vtkNew<vtkCellArray> triangles;
        for( auto c : cells ){
            auto infv = c->index(T.infinite_vertex());
            triangles->InsertNextCell(3);
            for( auto j=0; j < 4; ++j){
                if (j == infv)
                    continue;
                triangles->InsertCellPoint(c->vertex(j)->info());
            }
        }
        poly->SetPolys(triangles);
        vtkNew<vtkPolyDataWriter> writer;
        writer->SetFileName("Mesh.vtk");
        writer->SetInputData(poly);
        writer->Write();
        */
    }
    float diff((float)clock() - (float)t1);
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    return 0;
}
