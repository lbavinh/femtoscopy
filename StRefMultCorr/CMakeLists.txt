# CMakeLists.txt for StFemtoEvent package. It creates a library with dictionary
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(StRefMultCorr)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
if(ROOT_CMAKE)
  include(${ROOT_USE_FILE})
endif(ROOT_CMAKE)

#---Define _VANILLA_ROOT_ variable needed for the project
add_definitions(-D_VANILLA_ROOT_)

set(INCLUDE_DIRECTORIES
  ${CMAKE_SOURCE_DIR}
  ${ROOT_INCLUDE_DIRS}
)

include_directories(${INCLUDE_DIRECTORIES})
add_definitions(${ROOT_CXX_FLAGS})

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "-O2")

set(HEADER_FILES
	"${CMAKE_CURRENT_SOURCE_DIR}/CentralityMaker.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/StRefMultCorr.h"
	# "${CMAKE_CURRENT_SOURCE_DIR}/BadRun.h"
	"${CMAKE_CURRENT_SOURCE_DIR}/Param.h"
)

set(CXX_FILES
	"${CMAKE_CURRENT_SOURCE_DIR}/CentralityMaker.cxx"
	"${CMAKE_CURRENT_SOURCE_DIR}/StRefMultCorr.cxx"
	"${CMAKE_CURRENT_SOURCE_DIR}/Param.cxx"

)

#---Generate dictionary for the project
ROOT_GENERATE_DICTIONARY(G__StRefMultCorr
	${HEADER_FILES}
	LINKDEF StRefMultCorrLinkDef.h)

#---Create a shared library with generated dictionary
add_library(StRefMultCorr SHARED
	${CXX_FILES}
	G__StRefMultCorr.cxx)

target_link_libraries(StRefMultCorr ${ROOT_LIBRARIES})
