cmake_minimum_required(VERSION 2.8)
project(proto-sph-1)

if ( NOT WIN32 )
    # does not work for MSVC...?
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
    if(COMPILER_SUPPORTS_CXX11)
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has C++11 support.")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    elseif(COMPILER_SUPPORTS_CXX0X)
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has C++0x support.")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    else()
        message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
    endif()
endif ()

if ( WIN32 )
    set (CMAKE_PREFIX_PATH "C:\\Qt\\5.5\\msvc2013_64\\")
endif ()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../sceneparser/src")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../particlegenerator/src")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../particlegenerator" "${CMAKE_CURRENT_BINARY_DIR}/particlegenerator")

add_definitions(-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

find_package(OpenMP)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

find_package(Boost REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

find_library(H5Part H5Part REQUIRED)
find_library(HDF5 hdf5 REQUIRED)
find_library(QGLVIEWER QGLViewer REQUIRED)

option(NOGUI OFF)

if (${NOGUI})
    find_package(Qt5Core REQUIRED)
    add_definitions(-DNOGUI)
    include_directories(${Qt5Core_INCLUDE_DIRS})

    # Add compiler flags for building executables (-fPIE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Core_EXECUTABLE_COMPILE_FLAGS}")
else ((${NOGUI}))
    find_package(Qt5Core REQUIRED)
    include_directories(${Qt5Core_INCLUDE_DIRS})
    find_package(Qt5Gui REQUIRED)
    include_directories(${Qt5Gui_INCLUDE_DIRS})
    find_package(Qt5Xml REQUIRED)
    include_directories(${Qt5Xml_INCLUDE_DIRS})
    find_package(Qt5PrintSupport REQUIRED)
    include_directories(${Qt5PrintSupport_INCLUDE_DIRS})
    find_package(Qt5OpenGL)
    include_directories(${Qt5OpenGL_INCLUDE_DIRS})
    set(QT_USE_QTOPENGL 1)
    find_package(OpenGL REQUIRED)
    include_directories(${OpenGL_INCLUDE_DIRS})

    # Add compiler flags for building executables (-fPIE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Core_EXECUTABLE_COMPILE_FLAGS} ${Qt5Gui_EXECUTABLE_COMPILE_FLAGS} ${Qt5OpenGL_EXECUTABLE_COMPILE_FLAGS}")
endif (${NOGUI})

add_subdirectory(src)

#testing stuff

option(WITH_TESTING OFF)

if (${WITH_TESTING})
    add_subdirectory(tests)
endif (${WITH_TESTING})

message(STATUS ${CMAKE_CXX_COMPILER_ID})
