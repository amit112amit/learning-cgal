#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <Eigen/Dense>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned,K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef Tds::Vertex_handle Vertex_handle;
typedef Tds::Cell_handle Cell_handle;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef Delaunay::Point Point;
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

        rotMat = cos_t*Matrix3d::Identity() + sin_t*axis_cross +
                (1-cos_t)*outer;
        Matrix3Xd rPts(3,N);
        rPts = rotMat*points; // The points on a sphere rotated
	points = rPts;

        std::vector<std::pair<Point,unsigned>> spherePoints;
        spherePoints.push_back(std::make_pair(Point(0.,0.,0.),N));
        for( auto i=0; i < N; ++i){
            spherePoints.push_back( std::make_pair(Point(points(0,i),
                                                   points(1,i),
                                                   points(2,i)),
                                             i));
        }

        // Calculate the convex hull
        Delaunay T(spherePoints.begin(),spherePoints.end());

        // To extract the surface
        std::vector<Cell_handle> cells;
        T.tds().incident_cells( T.infinite_vertex(),
                                std::back_inserter(cells) );

        // Write to a vtk file
        vtkNew<vtkCellArray> triangles;
        for( auto c : cells ){
            auto infv = c->index(T.infinite_vertex());
            //triangles->InsertNextCell(3);
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
    }
    float diff(static_cast<float>(clock()) - static_cast<float>(t1));
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    return 0;
}
