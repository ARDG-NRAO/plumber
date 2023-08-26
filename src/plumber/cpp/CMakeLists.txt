project(ZernikePolynomials LANGUAGES CXX)

cmake_minimum_required(VERSION 3.12.0)
set (CMAKE_CXX_STANDARD 11)
set(OPT "-g -O2 ")
include(FetchContent)

set(CMAKE_C_FLAGS "${OPT} -std=c++11 -lm -Wno-deprecated \
           ${DEFS64}                                        \
           -Wall -pipe -Wno-non-template-friend             \
           -Woverloaded-virtual -Wcast-align                \
           -Wno-comment -D_FILE_OFFSET_BITS=64              \
           -D_LARGEFILE_SOURCE ")

set(INSTALLDIR "install")

FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
        FetchContent_Populate(
        pybind11
        URL https://github.com/pybind/pybind11/archive/refs/tags/v2.10.4.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pybind11
        )          
        add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pybind11)        
endif()
set(WRAPPER_LIST ZernikePolyCalc2py)

foreach(WRAPPER IN LISTS WRAPPER_LIST)
    pybind11_add_module(${WRAPPER} ${WRAPPER}.cc)
    target_link_libraries(${WRAPPER} PRIVATE ${APP_LINK_LIBRARIES})
    install(TARGETS ${WRAPPER} DESTINATION "${CMAKE_SOURCE_DIR}/install") 
endforeach()
