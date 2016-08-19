cmake_minimum_required(VERSION 3.5)
project(energy_alpscore)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(ALPSCore COMPONENTS hdf5 accumulators mc params gf REQUIRED)
message(STATUS "ALPSCore includes: " ${ALPSCore_INCLUDE_DIR})
include_directories(${ALPSCore_INCLUDE_DIR})

find_package(TCLAP)
message(STATUS "TCLAP includes: " ${TCLAP_INCLUDE_DIR})
include_directories(${TCLAP_INCLUDE_DIR})

find_package(GSL)
message(STATUS "GSL includes: " ${GSL_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(energy_alpscore ${SOURCE_FILES})
add_executable(tail_check tail_check.cpp)

target_link_libraries(energy_alpscore
       ${ALPSCore_LIBRARIES} ${GSL_LIBRARIES}
        )

target_link_libraries(tail_check
        ${ALPSCore_LIBRARIES} ${GSL_LIBRARIES}
        )