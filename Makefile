######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=g++
GCC=gcc

GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
#CFLAGS=-std=c++0x -Wall -DVERSION=\"$(GIT_VERSION)\" 
CFLAGS= -g -Wall -DVERSION=\"$(GIT_VERSION)\" -std=c++0x #-D_NO_RAZF
INCLUDE=-Isrc/lib -Isrc/bamtools/include -Isrc/bamtools/src -Isrc/Complete-Striped-Smith-Waterman-Library/src/ -Isrc/seqan/core/include/ -Isrc/seqan/extras/include 
OUTFOLD=bin/
LIBS=-L./ -lbamtools -fopenmp -lz -lm
RUNTIME=-Wl,-rpath=src/bamtools/lib/

all: createBin bamtools libbamtools.a buildWHAMBAM clean
debug: createBin bamtools libbamtools.a buildWHAMBAMD clean

createBin:
	-mkdir bin
bamtools:
	cd src/bamtools && mkdir -p build && cd build && cmake .. && make
libbamtools.a: bamtools
	cp src/bamtools/lib/libbamtools.a .
razf.o:
	cd src/lib && $(GCC) -c razf.c && cd ../../

cpsw:
	cp src/Complete-Striped-Smith-Waterman-Library/src/*_cpp* src/lib

buildWHAMBAM: libbamtools.a 
	$(CC) $(CFLAGS) src/lib/*cpp  src/Complete-Striped-Smith-Waterman-Library/src/*cpp src/fastahack/*cpp src/bin/multi-wham-testing.cpp $(INCLUDE) $(LIBS)  -o $(OUTFOLD)WHAM-BAM $(RUNTIME)
buildWHAMBAMD: libbamtools.a razf.o
	$(CC) $(CFLAGS) -g -DDEBUG src/lib/*cpp src/bin/multi-wham-testing.cpp $(INCLUDE) $(LIBS)  -o $(OUTFOLD)WHAM-BAM $(RUNTIME)
buildWHAMDUMPER:
	$(CC) $(CFLAGS) -g src/lib/*cpp   src/bin/multi-wham.cpp $(INCLUDE) $(LIBS) -o $(OUTFOLD)WHAM-BAM-DUMPER $(RUNTIME)
buildWHAMBAMGENE:
	$(CC) $(CFLAGS) -g src/lib/*cpp  src/bin/multi-wham-testing-gene.cpp  $(INCLUDE) $(LIBS) -o $(OUTFOLD)WHAM-BAM-GENE $(RUNTIME)

clean:
	-@rm *.a
