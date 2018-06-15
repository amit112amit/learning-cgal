# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.
cmake_minimum_required(VERSION 2.8.11)
project( CPP CXX )

# CGAL and its components
find_package( CGAL REQUIRED )

# include helper file
include( ${CGAL_USE_FILE} )

# CGAL and its components
find_package( VTK REQUIRED )

# include helper file
include( ${VTK_USE_FILE} )

# Boost and its components
find_package( Boost REQUIRED )

# Eigen
find_package( Eigen3 3.3 REQUIRED NO_MODULE )

add_executable(cgalStereo cgalStereo.cxx)
add_executable(cgalData cgalData.cxx)
add_executable(cgal3d cgal3d.cxx)
add_executable(cgal3dCH cgal3dCH.cxx)
add_executable(vtk3d vtk3d.cxx)
add_executable(vtk2d vtkStereo.cxx)
add_executable(simple simple.cxx)

set( programs cgalStereo cgal3d cgal3dCH vtk3d vtk2d cgalData )

foreach(program ${programs})
        target_link_libraries( ${program} PRIVATE Eigen3::Eigen ${VTK_LIBRARIES} ${CGAL_LIBRARIES})
        target_include_directories( ${program} PRIVATE ${VTK_INCLUDE_DIRS} ${CGAL_INCLUDE_DIRS} )
endforeach(program)
