# This is the cmake list for all wham tools 
# Author zev.kronenberg@gmail.com
# 2013

cmake_minimum_required(VERSION 2.6)	

project(WHAM CXX)

set(CMAKE_CXX_FLAGS "-lboost_program_options")
#set(CMAKE_CXX_FLAGS "-std=c++0x -lboost_program_options")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(SRC ${CMAKE_SOURCE_DIR}/src/lib/flag.cpp ${CMAKE_SOURCE_DIR}/src/lib/read_pileup.cpp)
set(HEADERS ${CMAKE_SOURCE_DIR}/src/lib/flag.h ${CMAKE_SOURCE_DIR}/src/lib/read_pileup.h)

find_package(Boost  REQUIRED) 
include_directories(${Boost_INCLUDE_DIRS}) 


add_library(lib ${SRC} ${HEADERS})

include_directories(/usr/local/include/bamtools)
include_directories(${CMAKE_SOURCE_DIR}/src/lib)

link_directories(/usr/local/lib/bamtools/)
link_directories(${CMAKE_SOURCE_DIR}/src/lib)

add_executable(wham ${CMAKE_SOURCE_DIR}/src/bin/wham.cpp  ${SRC} ${HEADER})
add_executable(whamy ${CMAKE_SOURCE_DIR}/src/bin/multi-wham.cpp  ${SRC} ${HEADER})
add_executable(emissions ${CMAKE_SOURCE_DIR}/src/bin/hmm_region.cpp  ${SRC} ${HEADER})
add_executable(raw ${CMAKE_SOURCE_DIR}/src/bin/multi-wham-testing.cpp   ${SRC} ${HEADER})


target_link_libraries(wham ${lib}  bamtools)
target_link_libraries(whamy ${lib}  bamtools)
target_link_libraries(emissions ${lib}  bamtools)
target_link_libraries(raw ${lib} bamtools)







