cmake_minimum_required(VERSION 2.8.9)

include_directories(./../../include/simple)
file(GLOB SOURCES "*.cpp")

set(INCLUDE_DIRS ${GSL_INCLUDE_DIRS})

set(LIBS ${LIBS} ${GSL_LIBRARIES})

add_executable(siapop-simple ${SOURCES})

target_link_libraries(siapop-simple ${LIBS})
