######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=gcc
CXX=g++
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
CFLAGS= -fstack-protector-all -Wall -DVERSION=\"$(GIT_VERSION)\" -std=c++0x -Wno-sign-compare -O3
INCLUDE=-Isrc/lib -Isrc/bamtools/include -Isrc/bamtools/src -Isrc/ -Isrc/fastahack -Isrc/Complete-Striped-Smith-Waterman-Library/src/
INCLUDE_PLUS := $(INCLUDE) -Isrc/seqan/core/include/ -Isrc/seqan/extras/include -Isrc/Complete-Striped-Smith-Waterman-Library/src/

SSW=src/Complete-Striped-Smith-Waterman-Library/src

LIBS=-fopenmp -lz -lm
LBAMTOOLS=src/bamtools/lib/libbamtools.a
CPP_FILES := $(wildcard src/lib/*.cpp)
OBJ_FILES := $(addprefix src/obj/,$(notdir $(CPP_FILES:.cpp=.o)))
OBJ_FILES := $(OBJ_FILES) src/fastahack/Fasta.o $(SSW)/ssw_cpp.o $(SSW)/ssw.o

all: createBin bin/whamg bin/wham

createBin:
	-mkdir bin
src/bamtools/lib/libbamtools.a:
	cd src/bamtools && mkdir -p build && cd build && cmake .. && make

src/Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.o:
	cd src/Complete-Striped-Smith-Waterman-Library/src && make

src/fastahack/Fasta.o:
	cd src/fastahack && make

src/obj/%.o: src/lib/%.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $(LIBS) -c -o $@ $<

bin/whamg: src/bamtools/lib/libbamtools.a $(OBJ_FILES)
	$(CXX) $(CFLAGS) $(INCLUDE) src/bin/whamg.cpp $(LBAMTOOLS) $(OBJ_FILES) -o bin/whamg $(LIBS)

bin/wham: src/bamtools/lib/libbamtools.a $(OBJ_FILES)
	$(CXX) $(CFLAGS) $(INCLUDE_PLUS) src/bin/wham.cpp $(LBAMTOOLS) $(OBJ_FILES) -o bin/wham $(LIBS)
clean:
	rm -rf bin && rm -rf src/obj/*