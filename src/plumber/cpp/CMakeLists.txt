project(ZernikePolynomials LANGUAGES CXX)

cmake_minimum_required(VERSION 3.24.0)
set (CMAKE_CXX_STANDARD 14)
set(OPT "-g -O2 ")
include(FetchContent)

set(CMAKE_C_FLAGS "${OPT} -std=c++14 -lm -Wno-deprecated \
           ${DEFS64}                                        \
           -Wall -pipe -Wno-non-template-friend             \
           -Woverloaded-virtual -Wcast-align                \
           -Wno-comment -D_FILE_OFFSET_BITS=64              \
           -D_LARGEFILE_SOURCE ")

set(INSTALLDIR "install")
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/${INSTALLDIR})
set(CMAKE_INSTALL_RPATH ${CMAKE_SOURCE_DIR}/${INSTALLDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(LIBDIR ${CMAKE_SOURCE_DIR}/${INSTALLDIR}/lib)
set(INCDIR ${CMAKE_SOURCE_DIR}/${INSTALLDIR}/include)

FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.11.1
    GIT_SHALLOW FALSE)
FetchContent_MakeAvailable(pybind11)

FetchContent_Declare(xtl
  GIT_REPOSITORY https://github.com/xtensor-stack/xtl
  GIT_TAG 0.7.5
  GIT_SHALLOW FALSE)

FetchContent_MakeAvailable(xtl)

FetchContent_Declare(xtensor
  GIT_REPOSITORY https://github.com/xtensor-stack/xtensor
  GIT_TAG 0.24.6
  GIT_SHALLOW FALSE)

FetchContent_MakeAvailable(xtensor)

FetchContent_Declare(xtensor-python
    GIT_REPOSITORY https://github.com/xtensor-stack/xtensor-python
    GIT_TAG 0.24.1
    GIT_SHALLOW FALSE)
FetchContent_MakeAvailable(xtensor-python)

set(WRAPPER_LIST ZernikePolyCalc2py)
add_library(ZernikePolyCalc SHARED ZernikePolyCalc.cc ZernikePolyCalc_xt.cc)
target_link_libraries(ZernikePolyCalc PRIVATE xtensor)

foreach(WRAPPER IN LISTS WRAPPER_LIST)
    pybind11_add_module(${WRAPPER} ${WRAPPER}.cc)
    # target_link_libraries(${WRAPPER} PRIVATE ${APP_LINK_LIBRARIES})
    install(TARGETS ${WRAPPER} DESTINATION "${CMAKE_SOURCE_DIR}/install") 
endforeach()

install(TARGETS ZernikePolyCalc DESTINATION "${CMAKE_SOURCE_DIR}/install") 
