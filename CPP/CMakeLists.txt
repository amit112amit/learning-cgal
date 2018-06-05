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

# Creating entries for all C++ files with "main" routine
# ##########################################################
include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program("cgalStereo.cxx")

# Eigen

find_package( Eigen3 3.3 REQUIRED NO_MODULE )

target_link_libraries( cgalStereo PRIVATE Eigen3::Eigen ${VTK_LIBRARIES})
target_include_directories( cgalStereo PRIVATE ${VTK_INCLUDE_DIRS} )

add_executable(vtk3d vtk3d.cxx)
target_link_libraries( vtk3d PRIVATE Eigen3::Eigen ${VTK_LIBRARIES})
target_include_directories( vtk3d PRIVATE ${VTK_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS})

add_executable(vtk2d vtkStereo.cxx)
target_link_libraries( vtk2d PRIVATE Eigen3::Eigen ${VTK_LIBRARIES})
target_include_directories( vtk2d PRIVATE ${VTK_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS})
