cmake_minimum_required(VERSION 2.8.9)

include_directories(./../../include/time-dependent)
file(GLOB SOURCES "*.cpp")

set(INCLUDE_DIRS ${GSL_INCLUDE_DIRS})

set(LIBS ${LIBS} ${GSL_LIBRARIES})

add_executable(siapop-td ${SOURCES})

target_link_libraries(siapop-td ${LIBS})
