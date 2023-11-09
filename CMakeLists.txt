# Config file for compilation with CMake
# Only relevant for Batch mode, but this file needs to live in the RScore folder

add_library(RScore Species.cpp Cell.cpp Community.cpp FractalGenerator.cpp Genome.cpp Individual.cpp Landscape.cpp Model.cpp Parameters.cpp Patch.cpp Population.cpp RandomCheck.cpp RSrandom.cpp SubCommunity.cpp)

# pass config definitions to compiler
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	target_compile_definitions(RScore PRIVATE "RSDEBUG" "RSWIN64" "LINUX_CLUSTER")
else() # Windows
	target_compile_definitions(RScore PRIVATE "RSDEBUG" "RSWIN64")
endif()