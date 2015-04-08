cmake_minimum_required(VERSION 2.8)
set(ROS_BUILD_TYPE Release)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)
rosbuild_init()
include_directories(${PROJECT_SOURCE_DIR}/code)
#PROJECT( sbpl )
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-std=c++0x")
# Add a ROS C++ macro to tell SBPL to use ROS print statements
add_definitions(-DROS)

find_package( PkgConfig REQUIRED)
pkg_check_modules( gsl REQUIRED gsl )

rosbuild_add_library(sbpl
  src/discrete_space_information/environment_nav2D.cpp
  src/discrete_space_information/environment_navxythetalat.cpp
  src/discrete_space_information/environment_navxythetamlevlat.cpp
  src/discrete_space_information/environment_nav2Duu.cpp
  src/discrete_space_information/environment_XXX.cpp
  src/discrete_space_information/environment_robarm.cpp
  src/planners/araplanner.cpp
  src/planners/adplanner.cpp
  src/planners/viplanner.cpp
  src/planners/ppcpplanner.cpp
  src/planners/ANAplanner.cpp
  src/planners/rstarplanner.cpp
  src/planners/mha_planner.cpp
  src/planners/lazyARA.cpp
  src/utils/heap.cpp
  src/utils/mdp.cpp
  src/utils/utils.cpp
  src/utils/2Dgridsearch.cpp
  src/utils/config.cpp	
  )
target_link_libraries(sbpl ${gsl_LIBRARIES})

#set(SBPL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
#set(SBPL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
  
# Build executables in the bin directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

rosbuild_add_executable(test_sbpl src/test/main.cpp)
target_link_libraries(test_sbpl sbpl)

rosbuild_add_executable(test_adjacency_list src/test/test_adjacency_list.cpp)
target_link_libraries(test_adjacency_list sbpl)

# Download needed data file
rosbuild_download_test_data(http://pr.willowgarage.com/data/sbpl/env_examples/nav2d/willow-25mm-inflated-env.cfg env_examples/nav2d/willow-25mm-inflated-env.cfg 1b4551eeba6df03b6e597e598bbaa4f3)
