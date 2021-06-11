cmake_minimum_required(VERSION 3.11)

project(render
  VERSION 0.0.1
  LANGUAGES CXX)

aux_source_directory(src CXX_SOURCES)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

find_package(OpenMP REQUIRED)

add_executable(render "${CXX_SOURCES}")

target_link_libraries(render PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(render PUBLIC src)

