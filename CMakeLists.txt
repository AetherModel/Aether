cmake_minimum_required(VERSION 3.10)

project(Aether LANGUAGES CXX
               VERSION   2021.0)

# Directory variables
file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
set(MAIN_DIR ${PROJECT_SOURCE_DIR}/src/main)
set(TESTS_DIR ${PROJECT_SOURCE_DIR}/tests)
set(OUT_DIR ${PROJECT_SOURCE_DIR}/src/output)
set(RUN_DIR ${PROJECT_SOURCE_DIR}/run)

# Executable decision
if(USE_NETCDF)
  set(SRC_FILES ${OUT_DIR}/output_netcdf.cpp ${SRC_FILES})
else()
  set(SRC_FILES ${OUT_DIR}/output_binary.cpp ${SRC_FILES})
endif()
add_executable(aether ${SRC_FILES} ${MAIN_DIR}/main.cpp)
if(USE_DOUBLE_PRECISION)
  add_compile_definitions(AETHER_USE_PRECISION_DOUBLE)
endif()

# Tests executable
add_executable(tests ${TESTS_DIR}/electrodynamics_unit_tests.cpp)

# Set compiler and cmake options
target_compile_options(aether PUBLIC -ffast-math -O3)
set(SHARE_INCLUDE ${PROJECT_SOURCE_DIR}/share/include)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_FIND_ROOT_PATH ${SHARE_INCLUDE}
                         ${CMAKE_FIND_ROOT_PATH}
                         )
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}
                      ${PROJECT_SOURCE_DIR}/share
                      )
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/share
                      ${CMAKE_MODULE_PATH})
target_include_directories(aether PUBLIC ${PROJECT_SOURCE_DIR}/include)
include_directories(${SHARE_INCLUDE})

# Find dependencies
find_package(Armadillo REQUIRED)
target_include_directories(aether PUBLIC ${ARMADILLO_INCLUDE_DIRS})
find_package(BLAS)
find_package(LAPACK)
if(BLAS_FOUND AND LAPACK_FOUND)
  set(ARMA_DONT_USE_WRAPPER True)
  set(ARMA_USE_BLAS True)
  set(ARMA_USE_LAPACK True)
  target_link_libraries(aether PUBLIC ${BLAS_LIBRARIES})
  target_link_libraries(aether PUBLIC ${LAPACK_LIBRARIES})
elseif(ARMADILLO_FOUND)
  message(WARNING "BLAS and LAPACK not found. Matrix operation will be slower.")
  target_link_libraries(aether PUBLIC ${ARMADILLO_LIBRARIES})
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
  target_link_libraries(aether PUBLIC OpenMP::OpenMP_CXX)
endif()

if(USE_NETCDF)
  set(NETCDF_CXX "YES")
  find_package(NetCDF REQUIRED)
  target_include_directories(aether PUBLIC ${NETCDF_INCLUDES})
  target_link_libraries(aether PUBLIC ${NETCDF_LIBRARIES_CXX})
  target_link_libraries(aether PUBLIC ${NETCDF_LIBRARIES_C})
endif()

# Set up links for executable
file(CREATE_LINK ${CMAKE_BINARY_DIR}/aether
                 ${PROJECT_SOURCE_DIR}/share/run/aether SYMBOLIC)

