cmake_minimum_required(VERSION 2.8.3)
project(sbpl)
find_package(catkin REQUIRED COMPONENTS roscpp sbpl_utils)

set(ROS_BUILD_TYPE Release)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-std=c++0x")
# add_definitions(-DROS)

find_package( PkgConfig REQUIRED)
pkg_check_modules( gsl REQUIRED gsl )

include_directories(${PROJECT_SOURCE_DIR}/include ${catkin_INCLUDE_DIRS})

catkin_package(
  CATKIN_DEPENDS roscpp sbpl_utils
  DEPENDS gsl
  INCLUDE_DIRS include
  LIBRARIES sbpl_ros
  )

add_library(sbpl_ros
  src/discrete_space_information/environment_nav2D.cpp
  src/discrete_space_information/environment_navxythetalat_mha.cpp
  # src/discrete_space_information/environment_navxythetamlevlat.cpp
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
target_link_libraries(sbpl_ros ${gsl_LIBRARIES} ${catkin_LIBRARIES})

# Download needed data file
download_test_data(http://pr.willowgarage.com/data/sbpl/env_examples/nav2d/willow-25mm-inflated-env.cfg 
  env_examples/nav2d/willow-25mm-inflated-env.cfg 1b4551eeba6df03b6e597e598bbaa4f3)

#set(SBPL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
#set(SBPL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")

# add_executable(test_sbpl src/test/main.cpp)
# target_link_libraries(test_sbpl sbpl_ros)

add_executable(test_pa_islands src/test/main_pa_islands.cpp)
target_link_libraries(test_pa_islands sbpl_ros)

add_executable(test_adjacency_list src/test/test_adjacency_list.cpp)
target_link_libraries(test_adjacency_list sbpl_ros)
