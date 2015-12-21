cmake_minimum_required(VERSION 2.8)
project(particlegenerator)

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

find_package(Qt5Core REQUIRED)
include_directories(${Qt5Core_INCLUDE_DIRS})
add_definitions(${Qt5Core_DEFINITIONS})
# Add compiler flags for building executables (-fPIE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Core_EXECUTABLE_COMPILE_FLAGS}")


include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../sceneparser/src")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../sceneparser" "${CMAKE_CURRENT_BINARY_DIR}/sceneparser")
add_subdirectory(src)

message(STATUS ${CMAKE_CXX_COMPILER_ID})
