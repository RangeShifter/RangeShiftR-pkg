# Config file for compilation with CMake

if(NOT batchmode) # that is, RScore as a standalone
	cmake_minimum_required(VERSION 3.10)
	# set the project name and version
	project(RScore VERSION 2.1.0)
	# specify the C++ standard

	set(CMAKE_CXX_STANDARD 20)
	set(CMAKE_CXX_STANDARD_REQUIRED True)

	add_executable(RScore Main.cpp Species.cpp Cell.cpp Community.cpp FractalGenerator.cpp GeneticFitnessTrait.cpp Individual.cpp Landscape.cpp Management.cpp Model.cpp NeutralStatsManager.cpp Parameters.cpp Patch.cpp Population.cpp DispersalTrait.cpp RSrandom.cpp NeutralTrait.cpp SpeciesTrait.cpp SubCommunity.cpp Utils.cpp "unit_tests/testIndividual.cpp" "unit_tests/testNeutralStats.cpp" "unit_tests/testPopulation.cpp")
else() # that is, RScore compiled as library within RangeShifter_batch

	add_library(RScore Species.cpp Cell.cpp Community.cpp FractalGenerator.cpp GeneticFitnessTrait.cpp Individual.cpp Landscape.cpp Management.cpp Model.cpp NeutralStatsManager.cpp Parameters.cpp Patch.cpp Population.cpp DispersalTrait.cpp RSrandom.cpp NeutralTrait.cpp SpeciesTrait.cpp SubCommunity.cpp Utils.cpp)
endif()

if(OMP)
	find_package(OpenMP COMPONENTS CXX)
	if(OpenMP_CXX_FOUND)
		target_link_libraries(RScore PUBLIC OpenMP::OpenMP_CXX)
	endif()
endif()

# enable LINUX_CLUSTER macro on Linux + macOS
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	add_compile_definitions("LINUX_CLUSTER")
endif()


if(NOT batchmode)
	target_include_directories(RScore PUBLIC "${PROJECT_BINARY_DIR}")
endif()
