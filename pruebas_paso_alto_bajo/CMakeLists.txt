cmake_minimum_required(VERSION 3.15)
project(KFR_Project)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -O2 -march=native)

# IMPORTANTE: activar multiarch para evitar problemas con avx2
set(KFR_ENABLE_MULTIARCH ON CACHE BOOL "" FORCE)

# Opcional pero recomendable (igual que en el otro proyecto)
set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)

add_subdirectory(kfr)

add_executable(filtro_kfr main.cpp)

# CAMBIO CLAVE: enlazar ambas librerías
target_link_libraries(filtro_kfr PRIVATE kfr kfr_dsp)