
# CMakeLists.txt for FAS package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(STARFlowAnalysis)

# You need to tell CMake where to find the ROOT installation. This can be done
# in a number of ways:
#   - ROOT built with classic configure/make use the provided
#   $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix
#   for ROOT

# find_package(Git)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net MathMore)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

add_definitions(${ROOT_CXX_FLAGS})

add_subdirectory(StPicoDst)
add_subdirectory(StRefMultCorr)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall")

set(FLOW_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/StPicoDst
  ${CMAKE_CURRENT_SOURCE_DIR}/StRefMultCorr
  ${ROOT_INLCUDE_DIRS}
)

include_directories(${FLOW_INCLUDE_DIRECTORIES})

set(FLOW_INCLUDE_LIBRARIES
  StPicoDst
  StRefMultCorr
  ${ROOT_LIBRARIES}
)

add_executable(PicoDstFemtoscopy PicoDstFemtoscopy.cpp)
target_link_libraries(PicoDstFemtoscopy ${FLOW_INCLUDE_LIBRARIES})
