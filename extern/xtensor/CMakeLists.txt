cmake_minimum_required(VERSION 3.15)

project(
  xtensor-wheel
  LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(xtensor
  GIT_REPOSITORY https://github.com/xtensor-stack/xtensor
  GIT_TAG ${SKBUILD_PROJECT_VERSION}
  GIT_SHALLOW FALSE)

FetchContent_MakeAvailable(xtensor)
