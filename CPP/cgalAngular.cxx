#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/IO/Triangulation_off_ostream_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <Eigen/Dense>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIdFilter.h>
#include <vtkSmartPointer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned,K> Vb;
typedef K::Point_2 Point;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Delaunay;
typedef Delaunay::Vertex_circulator Vertex_circulator;
typedef Delaunay::Line_face_circulator Line_face_circulator;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::Map<VectorXd,0,Eigen::InnerStride<3>> MapVXd;
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

        // Now we will calculate the angles theta, phi for all the points
        MapVXd x(points.data(),N), y(&(points(1,0)),N), z(&(points(2,0)),N);
        VectorXd phi(N), theta(N);
        theta = y.binaryExpr(x, [](double_t b, double_t a){ return std::atan2(b,a);} );
        phi = ((x.array().square() + y.array().square()).matrix().cwiseSqrt()
               ).binaryExpr(z, [](double_t b, double_t a){ return std::atan2(b,a);} );

        // Now we need to create points for triangulation
        std::vector< std::pair<Point,unsigned> > idPoints;
        for(auto j=0; j < N; ++j){
            idPoints.push_back(std::make_pair(Point(theta(j),phi(j)),j));
        }

        Delaunay dt(idPoints.begin(),idPoints.end());

        // Write the triangulation to a VTK file
        vtkNew<vtkCellArray> triangles;
        for( auto ffi = dt.finite_faces_begin(); ffi != dt.finite_faces_end(); ++ffi){
            triangles->InsertNextCell(3);
            for(auto j=0; j < 3; ++j)
                triangles->InsertCellPoint(ffi->vertex(j)->info());
        }

        // Iterate over vertices incident on infinite_faces and add new shifted points
        Vertex_circulator vc = dt.incident_vertices(dt.infinite_vertex()), done(vc);
        if (vc != 0) {
            do{
                auto pt = vc->point();
                if( pt.x() > 0 ){
                    auto ptTheta = pt.x() - 2*M_PI;
                    auto vh = dt.insert( Point(ptTheta, pt.y()));
                    vh->info() = vc->info();
                }
            }while(++vc != done);
        }

        // Now find all triangles that intersect the line theta=-M_PI and add to the triangles
        Point p1(-M_PI,0.0), p2(-M_PI,M_PI);
        Line_face_circulator lfc = dt.line_walk(p1,p2), doneLfc(lfc);
        if (lfc != 0) {
            do{
                if(!dt.is_infinite(lfc)){
                    triangles->InsertNextCell(3);
                    for(auto j=0; j < 3; ++j)
                        triangles->InsertCellPoint(lfc->vertex(j)->info());
                }
            }while(++lfc != doneLfc);
        }
        poly->SetPolys(triangles);

        vtkNew<vtkPolyDataWriter> writer;
        writer->SetInputData(poly);
        writer->SetFileName("Final.vtk");
        writer->Write();
    }

    float diff((float)clock() - (float)t1);
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    return 0;
}
