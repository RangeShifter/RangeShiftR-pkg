# Config file for compilation with CMake
# Only relevant for Batch mode, but this file needs to live in the RScore folder

add_library(RScore Species.cpp Cell.cpp Community.cpp FractalGenerator.cpp Genome.cpp Individual.cpp Landscape.cpp Model.cpp Parameters.cpp Patch.cpp Population.cpp RandomCheck.cpp RSrandom.cpp SubCommunity.cpp Utils.cpp)

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