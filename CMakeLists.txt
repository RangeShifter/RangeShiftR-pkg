# Config file for compilation with CMake
# Only relevant for Batch mode, but this file needs to live in the RScore folder

add_library(RScore Species.cpp Cell.cpp Community.cpp FractalGenerator.cpp Genome.cpp Individual.cpp Landscape.cpp Model.cpp Parameters.cpp Patch.cpp Population.cpp RandomCheck.cpp RSrandom.cpp SubCommunity.cpp)

# pass config definitions to compiler
target_compile_definitions(RScore PRIVATE "RSDEBUG" "LINUX_CLUSTER" "RSWIN64" "BATCH")