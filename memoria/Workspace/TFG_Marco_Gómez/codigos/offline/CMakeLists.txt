cmake_minimum_required(VERSION 3.15)
project(BOProject)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O2 -march=native -DNDEBUG)

add_compile_definitions(USE_NLOPT)

set(LIMBO_DIR "" CACHE PATH "Ruta al directorio src de Limbo")

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem thread)
find_package(TBB REQUIRED)
find_package(NLopt REQUIRED)

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

if(LIMBO_DIR)
    include_directories(${LIMBO_DIR}/src)
endif()

add_executable(consts_calculator consts_calculator.cpp scaling.cpp)
target_link_libraries(consts_calculator PRIVATE kfr kfr_dsp)

find_package(yaml-cpp REQUIRED)

add_executable(bo main.cpp scaling.cpp)
target_link_libraries(bo PRIVATE 
    kfr 
    kfr_dsp 
    Eigen3::Eigen 
    Boost::system 
    Boost::filesystem 
    Boost::thread 
    nlopt
    yaml-cpp
)