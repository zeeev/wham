######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=gcc
CXX=g++
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
CFLAGS= -fstack-protector-all -Wall -DVERSION=\"$(GIT_VERSION)\" -std=c++0x  -Wno-sign-compare
INCLUDE=-Isrc/lib -Isrc/bamtools/include -Isrc/bamtools/src -Isrc/ -Isrc/fastahack -Isrc/Complete-Striped-Smith-Waterman-Library/src/

LIBS=-fopenmp -lz -lm
LBAMTOOLS=src/bamtools/lib/libbamtools.a
CPP_FILES := $(wildcard src/lib/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
OBJ_FILES := $(OBJ_FILES) src/fastahack/Fasta.o

all: createBin bin/whamg

createBin:
	-mkdir bin
src/bamtools/lib/libbamtools.a:
	cd src/bamtools && mkdir -p build && cd build && cmake .. && make

src/fastahack/Fasta.o:
	cd src/fastahack && make

obj/%.o: src/lib/%.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) $(LIBS) -c -o $@ $<

bin/whamg: src/bamtools/lib/libbamtools.a $(OBJ_FILES)
	$(CXX) $(CFLAGS) $(INCLUDE) src/bin/whamg.cpp $(LBAMTOOLS) $(OBJ_FILES) -o bin/whamg $(LIBS)

clean:
	rm -rf bin && rm -rf obj/*