cmake_minimum_required(VERSION 2.8.9)

# include header files
include_directories(./../../include/constant-rate)
# include all source files
file(GLOB SOURCES "*.cpp")

set(INCLUDE_DIRS ${GSL_INCLUDE_DIRS})

set(LIBS ${LIBS} ${GSL_LIBRARIES})

add_executable(siapop ${SOURCES})

target_link_libraries(siapop ${LIBS})
