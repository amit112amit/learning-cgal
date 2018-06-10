#include <Eigen/Dense>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIdFilter.h>
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
        vtkNew<vtkPolyDataReader> reader;
        reader->SetFileName("T7.vtk");
        reader->Update();
        auto poly = reader->GetOutput();
        auto N = poly->GetNumberOfPoints();
        auto pts = (double*) poly->GetPoints()->GetData()->GetVoidPointer(0);
        Map3Xd points(pts,3,N);

        // Project points to unit sphere
        points.colwise().normalize();

        // Create a 3D tessellation
        auto idf = vtkSmartPointer<vtkIdFilter>::New();
        idf->PointIdsOn();
        idf->SetIdsArrayName("OrigIds");
        idf->SetInputData(poly);
        auto d3d = vtkSmartPointer<vtkDelaunay3D>::New();
        d3d->SetInputConnection(idf->GetOutputPort());
        auto dssf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        dssf->SetInputConnection(d3d->GetOutputPort());
        dssf->Update();
        auto final = dssf->GetOutput();
        auto interim = final->GetPolys();
        interim->InitTraversal();
        auto origIds = vtkIdTypeArray::SafeDownCast(
                                final->GetPointData()->GetArray("OrigIds"));
        auto pointIds = vtkSmartPointer<vtkIdList>::New();
        auto finalCells = vtkSmartPointer<vtkCellArray>::New();
        while(interim->GetNextCell(pointIds)){
            int numIds = pointIds->GetNumberOfIds();
            finalCells->InsertNextCell(numIds);
            for(auto j=0; j < numIds; j++ ){
                int id = (int)origIds->GetTuple1( pointIds->GetId(j) );
                finalCells->InsertCellPoint(id);
            }
        }
        poly->SetPolys(finalCells);
        auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName("Mesh.vtk");
        writer->SetInputData(poly);
        writer->Write();

    }

    float diff((float)clock() - (float)t1);
    std::cout << "Time elapsed : " << diff / CLOCKS_PER_SEC
              << " seconds" << std::endl;
    return 0;
}
