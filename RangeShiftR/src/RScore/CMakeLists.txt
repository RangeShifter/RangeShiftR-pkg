# Config file for compilation with CMake

if (NOT batchmode) # that is, RScore as a standalone
	cmake_minimum_required(VERSION 3.10)
	# set the project name and version
	project(RScore VERSION 2.1.0)
	# specify the C++ standard
	set(CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_STANDARD_REQUIRED True)
	add_executable(RScore Main.cpp Species.cpp Cell.cpp Community.cpp FractalGenerator.cpp Genome.cpp Individual.cpp Landscape.cpp Model.cpp Parameters.cpp Patch.cpp Population.cpp RandomCheck.cpp RSrandom.cpp SubCommunity.cpp Utils.cpp)
else() # that is, RScore compiled as library within RangeShifter_batch
	add_library(RScore Species.cpp Cell.cpp Community.cpp FractalGenerator.cpp Genome.cpp Individual.cpp Landscape.cpp Model.cpp Parameters.cpp Patch.cpp Population.cpp RandomCheck.cpp RSrandom.cpp SubCommunity.cpp Utils.cpp)
endif()

# pass config definitions to compiler
target_compile_definitions(RScore PRIVATE RSWIN64)

# enable LINUX_CLUSTER macro on Linux + macOS
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	add_compile_definitions("LINUX_CLUSTER")
endif()

# Debug Mode by default, unless "release" is passed
if(NOT DEFINED release)
	add_compile_definitions(RSDEBUG)
endif()

if(NOT batchmode)
	target_include_directories(RScore PUBLIC "${PROJECT_BINARY_DIR}")
endif() 