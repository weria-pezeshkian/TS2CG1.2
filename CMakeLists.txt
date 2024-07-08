cmake_minimum_required(VERSION 3.10)
project(TS2CG)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(TS2CG/Solvate)
add_subdirectory(TS2CG/Pointillism)
add_subdirectory(TS2CG/MembraneBuilder)

install(TARGETS SOL PLM PCG
        RUNTIME DESTINATION TS2CG)
