cmake_minimum_required(VERSION 3.20)
include(FindPackageHandleStandardArgs)
# ---------------------------------------------------------------------
project(COULCC-25 VERSION 0.9
        DESCRIPTION ""
        LANGUAGES Fortran)
# ---------------------------------------------------------------------
enable_language(Fortran)
# ---------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)
set(CMAKE_INSTALL_PREFIX ".") 
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/mod)
# ---------------------------------------------------------------------
add_library(COULCC  coulcc.f90)
add_executable(testcc test/test-program.f90)
target_link_libraries(testcc PUBLIC COULCC)
