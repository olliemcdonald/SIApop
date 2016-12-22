cmake_minimum_required(VERSION 2.8.9)
project(siapop)

# Set verbose for generating makefile and c++ flags
set(CMAKE_VERBOSE_MAKEFILE ON)
set(MAKE_CXX_FLAGS "-g -Wall")

# executables placed in ./bin and libraries in ./lib
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/lib)

# Include macro to find gsl library location
find_package(GSL REQUIRED)
if (NOT GSL_FOUND)
  message(FATAL_ERROR "Package GSL not found")
endif (NOT GSL_FOUND)

set(INCLUDE_DIRS ${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

# add subdirectories that contain own CMakeLists.txt
add_subdirectory(./src/constant-rate)
add_subdirectory(./src/simple)
add_subdirectory(./src/time-dependent)
