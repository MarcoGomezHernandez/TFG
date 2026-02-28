cmake_minimum_required(VERSION 3.15)
project(GeneticProject)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -O2 -march=native)

include(FetchContent)
FetchContent_Declare(
    kfr
    GIT_REPOSITORY https://github.com/kfrlib/kfr.git
    GIT_TAG        main
)

set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(kfr)

include_directories(/usr/local/Neun/0.4.0)

add_executable(consts_calculator consts_calculator.cpp scaling.cpp)
target_link_libraries(consts_calculator PRIVATE kfr kfr_dsp)

add_executable(genetic main.cpp scaling.cpp)
target_link_libraries(genetic PRIVATE kfr kfr_dsp)