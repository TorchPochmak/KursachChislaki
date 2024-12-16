cmake_minimum_required(VERSION 3.21)
project(nm_mtrx_3
        VERSION 2.0
        DESCRIPTION "Matrix library"
        LANGUAGES CXX)
        
add_library (
    nm_mtrx_3
    src/Matrix.cpp)

target_include_directories (
    nm_mtrx_3
    PUBLIC
    ./include)

set_target_properties(
    nm_mtrx_3 PROPERTIES
    LANGUAGES CXX
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    VERSION 2.0
    DESCRIPTION "Matrix")