cmake_minimum_required(VERSION 3.3)
project(HDP_experiments)

find_package(GSL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    corpus.cpp
    corpus.h
    hdp.cpp
    hdp.h
    likelihood_computer.cpp
    likelihood_computer.h
    main.cpp
    running.cpp
    running.h
    state.cpp
    state.h
    utils.cpp
    utils.h)

add_executable(HDP_experiments ${SOURCE_FILES})

if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    target_link_libraries (HDP_experiments ${GSL_LIBRARIES})
endif (GSL_FOUND)