cmake_minimum_required(VERSION 3.15)
# Project Name
if(SCEMENT_PROJECT_NAME)
#
#
project(
  ${SCEMENT_PROJECT_NAME}
  VERSION ${SCEMENT_PROJECT_VERSION}
  LANGUAGES CXX)
#
else()
#
set(SCEMENT_PROJECT_NAME scement)
set(SCEMENT_PROJECT_VERSION 1.0)
project(
    ${SCEMENT_PROJECT_NAME}
    VERSION ${SCEMENT_PROJECT_VERSION}
    LANGUAGES CXX)
#
endif()
# additional cmake files
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11/tools/")
#
message(STATUS "SCEMENT PROJ NAME :  ${SCEMENT_PROJECT_NAME}")
message(STATUS "SCEMENT PROJ VER  :  ${SCEMENT_PROJECT_VERSION}")
message(STATUS "PROJECT NAME      :  ${PROJECT_NAME}")
message(STATUS "PROJECT VERSION   :  ${PROJECT_VERSION}")
message(STATUS "MODULE PATH       :  ${CMAKE_MODULE_PATH}")
#
set(CMAKE_VERBOSE_MAKEFILE ON)
set(BUILD_SHARED_LIBS ON)
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
#
# include directories
include_directories("${PROJECT_SOURCE_DIR}/include/")
#
# Compile
set(PROJECT_PYMOD _${PROJECT_NAME}cpp)
set(PROJECT_INSTALL_DIR ${PROJECT_NAME})
#
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ext/pybind11)
#
#
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG)
# Python Module
python_add_library(${PROJECT_PYMOD} MODULE src/mod_scement.cpp WITH_SOABI)
target_link_libraries(${PROJECT_PYMOD} PRIVATE pybind11::headers)
#pybind11_add_module(${PROJECT_PYMOD} src/mod_scement.cpp)

#
add_compile_definitions(SCEMENT_VERSION=${PROJECT_VERSION})
add_compile_definitions(SCEMENT_MODULE_NAME=${PROJECT_PYMOD})

# OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP version:  ${OpenMP_CXX_VERSION}")
    include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_PYMOD} PUBLIC OpenMP::OpenMP_CXX)
endif()

# External header-only library: Eigen
add_library(Eigen INTERFACE)
target_include_directories(Eigen INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/eigen)
target_link_libraries(${PROJECT_PYMOD} PUBLIC Eigen)

#
# Armadillo
find_package(Armadillo)
if(Armadillo_FOUND)
    message(STATUS "Armadillo version:  ${ARMADILLO_VERSION_STRING} (${ARMADILLO_VERSION_NAME})")
    add_compile_definitions(ARMA_DONT_USE_WRAPPER)
    add_compile_definitions(SCEMENT_ARMA)
    include_directories(${ARMADILLO_INCLUDE_DIRS}) 
    target_link_libraries(${PROJECT_PYMOD} PUBLIC ${ARMADILLO_LIBRARIES})
endif()
#
# BLAS, LAPACK and LAPACKE
# BLAS
find_package(BLAS)
if(BLAS_FOUND)
    add_compile_definitions(ARMA_USE_BLAS)
    add_compile_definitions(EIGEN_USE_BLAS)
    target_link_libraries(${PROJECT_PYMOD} PUBLIC BLAS::BLAS)
endif()
# LAPACK
find_package(LAPACK)
if(LAPACK_FOUND)
    add_compile_definitions(ARMA_USE_LAPACK)
    target_link_libraries(${PROJECT_PYMOD} PUBLIC LAPACK::LAPACK)
endif()
# LAPACKE
find_package(LAPACKE)
if(LAPACKE_FOUND)
    add_compile_definitions(EIGEN_USE_LAPACKE)
    target_link_libraries(${PROJECT_PYMOD} PUBLIC ${LAPACKE_LIBRARIES})
endif()


install(TARGETS ${PROJECT_PYMOD} DESTINATION ${PROJECT_INSTALL_DIR})
