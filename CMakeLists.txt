# This is the cmake file for wham 
# Author zev.kronenberg@gmail.com
# 2013

cmake_minimum_required(VERSION 2.6)	

project(WHAM CXX)


set(SRC ${CMAKE_SOURCE_DIR}/lib/flag.cpp ${CMAKE_SOURCE_DIR}/lib/read_pileup.cpp)
set(HEADERS ${CMAKE_SOURCE_DIR}/lib/flag.h ${CMAKE_SOURCE_DIR}/lib/read_pileup.h)


add_library(lib ${SRC} ${HEADERS})

include_directories(/Users/zev/tools/bamtools/include)
include_directories(${CMAKE_SOURCE_DIR}/lib)

link_directories(/Users/zev/tools/bamtools/lib)
link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable(wham ${CMAKE_SOURCE_DIR}/bin/wham.cpp ${SRC} ${HEADER})

find_library(BAM
	     NAMES bamtools
             HINTS "/Users/zev/tools/bamtools/lib"
)

target_link_libraries( wham ${lib} bamtools)
