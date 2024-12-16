cmake_minimum_required(VERSION 3.21)
project (prsr_tst)

add_executable (
    prsr_tst
    test.cpp)

target_link_libraries (
    prsr_tst
    PUBLIC
    prsr)

set_target_properties(
        prsr_tst PROPERTIES
        LANGUAGES CXX
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        VERSION 1.0
        DESCRIPTION "Just test")