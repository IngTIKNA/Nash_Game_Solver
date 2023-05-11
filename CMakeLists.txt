cmake_minimum_required(VERSION 3.7)
project(NASH)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O2")

#set(CMAKE_CXX_FLAGS "-g -Ofast -I/usr/include/python3.6m -L/usr/local/lib -lpython3 -Werror=unknown-pragmas")


include_directories(libs/Eigen)

add_executable(NASH
        NashGame.cpp
        NashGame.h
        main.cpp)

find_package(Python COMPONENTS Development)
target_include_directories(NASH PRIVATE ${Python_INCLUDE_DIRS})
target_link_libraries(NASH ${Python_LIBRARIES})

