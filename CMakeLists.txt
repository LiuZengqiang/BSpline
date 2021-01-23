cmake_minimum_required(VERSION 3.17)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
#set(CMAKE_CXX_FLAGS_DEBUG "-std=c++14 -O3 -fopenmp -pthread")

project(BSpline)

find_package(OpenGL REQUIRED)
if(OpenGL_FOUND)
    include_directories(${OPENGL_INCLUDE_DIR})
endif()

find_package(GLUT REQUIRED)
if(GLUT_FOUND)
    include_directories(${GLUT_INCLUDE_DIR})
endif()

set(THIRD_PARTY_LIBS
        ${OPENGL_LIBRARY}
        ${GLUT_LIBRARY}
        )

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(./src)

add_executable(BSpline main.cpp)
target_link_libraries(BSpline BSpline_lib)