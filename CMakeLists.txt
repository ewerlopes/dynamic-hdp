cmake_minimum_required(VERSION 3.3)
project(HDP_experiments)

find_package(GSL)
if(NOT OPENGL_FOUND)
    message(ERROR "GSL library not found!")
endif(NOT OPENGL_FOUND)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS} include)

set(SOURCE_FILES
    src/corpus.cpp
    src/hdp.cpp
    src/hdp_dynamic.cpp
    src/likelihood_computer.cpp
    src/main.cpp
    src/running.cpp
    src/state.cpp
    src/utils.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries (${PROJECT_NAME} ${GSL_LIBRARIES})

