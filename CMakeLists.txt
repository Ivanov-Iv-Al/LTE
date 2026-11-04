cmake_minimum_required(VERSION 3.10)
project(LTE_PCI_Detector C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
else()
    set(CMAKE_C_FLAGS_DEBUG "-g -Wall -Wextra -pedantic")
endif()

find_library(MATH_LIB m)
if(NOT MATH_LIB)
    message(FATAL_ERROR "Math library not found!")
endif()

add_executable(lte_pci_detector src/main.c)

target_link_libraries(lte_pci_detector ${MATH_LIB})

if(EXISTS "${CMAKE_SOURCE_DIR}/data/signal.pcm")
    file(COPY "${CMAKE_SOURCE_DIR}/data/signal.pcm"
         DESTINATION "${CMAKE_BINARY_DIR}")
    message(STATUS "Copied signal.pcm to build directory")
else()
    message(WARNING "signal.pcm not found in data/ directory")
endif()


message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C compiler: ${CMAKE_C_COMPILER}")
message(STATUS "Math library: ${MATH_LIB}")