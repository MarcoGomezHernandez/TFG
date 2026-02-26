cmake_minimum_required(VERSION 3.10)
project(KFR_Project)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Desactivar multiarch simplifica el binario y acelera la compilación 
# si solo vas a ejecutarlo en tu laptop (opcional)
set(KFR_ENABLE_MULTIARCH OFF CACHE BOOL "" FORCE)

add_subdirectory(kfr)

add_executable(filtro_kfr main.cpp)

# Enlaza específicamente kfr_dsp para tener las funciones de filtros
target_link_libraries(filtro_kfr kfr_dsp)