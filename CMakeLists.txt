cmake_minimum_required (VERSION 2.6)
set( CMAKE_CXX_COMPILER "/usr/bin/g++" )
set( CMAKE_C_COMPILER "/usr/bin/gcc" )

set(PROG_NAME "genetic_progr")

project (${PROG_NAME})


set(BOOST_ROOT "$ENV{HOME}/distr/boost_1_64_0/")
set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost 1.64.0 REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIR})



find_package( PythonLibs 2.7 REQUIRED )
find_package( PythonInterp 2.7 REQUIRED )
include_directories(${PYTHON_INCLUDE_DIRS})


set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
#FILE(GLOB SOURCE ${PROJECT_SOURCE_DIR}/src/*.cpp)
set(SOURCE ${SRC_DIR}/tpd.cpp ${SRC_DIR}/ttk.cpp ${SRC_DIR}/tel.cpp ${SRC_DIR}/ttree.cpp ${SRC_DIR}/ttrees.cpp)
set(SOURCE ${SOURCE} ${SRC_DIR}/${PROG_NAME}.cpp)

FILE(GLOB OWN_HEADERS ${PROJECT_SOURCE_DIR}/include/*.h)

set(3_INC_DIR "${PROJECT_SOURCE_DIR}/../../3rdparty/include")
set(HEADERS
    "${3_INC_DIR}/matplotlibcpp.h"
)

include_directories(${PROJECT_SOURCE_DIR}/include/)
include_directories(${3_INC_DIR}/)

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set (_VERSION_MAJOR 1)
set (_VERSION_MINOR 0)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/_config.h.in"
  "${PROJECT_BINARY_DIR}/_config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")


set( CMAKE_VERBOSE_MAKEFILE off )


set(CMAKE_CXX_FLAGS "-Wall -std=c++11 -std=c++1y -O0 -pg")
#set(CMAKE_CXX_FLAGS "-Wno-deprecated -std=c++11 -std=c++1y")
if(Boost_FOUND)
    add_executable(${PROG_NAME} ${SOURCE} ${OWN_HEADERS} ${HEADERS}) 
    target_link_libraries( ${PROG_NAME} ${PYTHON_LIBRARIES} ${Boost_LIBRARIES}  )
endif()

