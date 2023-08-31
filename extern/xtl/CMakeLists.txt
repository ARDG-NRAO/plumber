cmake_minimum_required(VERSION 3.15)

project(
  xtl-wheel
  LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(xtl
  GIT_REPOSITORY https://github.com/xtensor-stack/xtl
  GIT_TAG ${SKBUILD_PROJECT_VERSION}
  GIT_SHALLOW FALSE)

FetchContent_MakeAvailable(xtl)
