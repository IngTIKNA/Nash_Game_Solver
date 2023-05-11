cmake_minimum_required(VERSION 3.11)
project(test_nash)
Include(FetchContent)


FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)

find_package(PkgConfig REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

include_directories(../libs/Eigen)

set(SRC	
        ../NashGame.cpp
        ../NashGame.h	
	test_nash.cpp        
	${eigen}
)

add_executable(${PROJECT_NAME} ${SRC})

#include_directories(${BASE_DIR}/src)
#include_directories(${BASE_DIR}/include)


target_compile_features(test_nash PRIVATE cxx_std_17)  
if (WIN32)
  target_compile_options(test_nash PRIVATE /WX /W4)
else()
  target_compile_options(test_nash PRIVATE -Wall)  
endif()

target_link_libraries(test_nash PRIVATE Catch2::Catch2)






