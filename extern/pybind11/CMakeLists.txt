cmake_minimum_required(VERSION 3.15)

project(
  pybind11
  VERSION "2.11.1"
  LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(pybind11
      GIT_REPOSITORY https://github.com/pybind11/pybind11.git
      GIT_TAG ${SKBUILD_PROJECT_VERSION}
      GIT_SHALLOW False)

FetchContent_MakeAvailable(pybind11)
