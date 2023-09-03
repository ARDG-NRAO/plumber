# FILEPATH: CMakeLists.txt

cmake_minimum_required(VERSION 3.10)

project(Plumber)

# Set the C++ standard to C++11
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)

# Set the output directory for the executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# Set the output directory for the libraries
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# Add the include directory to the include path
include_directories(${CMAKE_SOURCE_DIR}/include)

find_package(Python COMPONENTS NumPy)
if(Python_NumPy_FOUND)
    message(STATUS "NumPy found")
    include_directories(${Python_NumPy_INCLUDE_DIRS})
    link_directories(${Python_NumPy_LIBRARY_DIRS})
    add_definitions(${Python_NumPy_DEFINITIONS})
else()
    message(STATUS "NumPy not found")
endif()

# Add the source directory to the include path
add_subdirectory(extern/pybind11)

add_subdirectory(extern/xtl)

add_subdirectory(extern/xtensor)

add_subdirectory(extern/xtensor-python)


# Add the source files to a list
set(SOURCES
    src/plumber/cpp/ZernikePolyCalc_xt.cc
    src/plumber/cpp/ZernikePolyCalc_xt.h
)

set(LIBRARIES
    pybind11::module
    xtensor-python
    xtensor
    xtl)

# Build the executable
add_library(ZernikePolyCalc_xt SHARED ${SOURCES})
# Link the libraries to the executable
target_link_libraries(ZernikePolyCalc_xt PUBLIC ${LIBRARIES})
pybind11_add_module(ZernikePolyCalc_xt2py src/plumber/cpp/ZernikePolyCalc_xt2py.cc)
target_link_libraries(ZernikePolyCalc_xt2py PRIVATE ${LIBRARIES} ZernikePolyCalc_xt)
# Install the executable and libraries
install(TARGETS ZernikePolyCalc_xt DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(FILES ${LIBRARIES} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

install(TARGETS ZernikePolyCalc_xt2py DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
# install(FILES ${LIBRARIES} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
# install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
