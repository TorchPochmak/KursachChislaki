cmake_minimum_required(VERSION 3.21)
project(prsr
        VERSION 2.0
        DESCRIPTION "Parser"
        LANGUAGES CXX)
        
add_library (
    prsr
    src/Parser.cpp)

target_include_directories (
    prsr
    PUBLIC
    ./include)

set_target_properties(
    prsr PROPERTIES
    LANGUAGES CXX
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
    VERSION 2.0
    DESCRIPTION "Parser")