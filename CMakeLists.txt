cmake_minimum_required(VERSION 3.18)

cmake_policy(SET CMP0167 NEW)

project(josephson LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Boost REQUIRED COMPONENTS program_options filesystem system)
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

add_library(josephson_lib
    src/Config.cpp
    src/EnergyGrid.cpp
    src/AndreevAmplitude.cpp
    src/ScatteringMatrix.cpp
    src/JosephsonCurrent.cpp
    src/misc.cpp
)

target_include_directories(josephson_lib
    PUBLIC src
)

target_link_libraries(josephson_lib
    Eigen3::Eigen
    PUBLIC OpenMP::OpenMP_CXX
)

add_executable(josephson
    src/josephson.cpp
)

target_link_libraries(josephson
    PRIVATE
        josephson_lib
        Boost::program_options
)

pybind11_add_module(josephson_python
    src/bindings_python.cpp
)

target_link_libraries(josephson_python
    PRIVATE josephson_lib
)

set_target_properties(josephson_python PROPERTIES
    OUTPUT_NAME josephson
)

add_executable(josephson_test
    src/josephson_test.cpp
)

target_link_libraries(josephson_test
    PRIVATE
        josephson_lib
        Eigen3::Eigen
        Boost::filesystem
        Boost::system
        Boost::program_options
)

enable_testing()

add_test(
    NAME josephson_test
    COMMAND josephson_test
)

add_test(
    NAME python_tests
    COMMAND ${CMAKE_COMMAND} -E env
        PYTHONPATH=${CMAKE_BINARY_DIR}
        pytest ${CMAKE_SOURCE_DIR}/python/test
)