######################################

# Makefile written by Zev Kronenberg #

#     zev.kronenberg@gmail.com       #

######################################



CC=g++
GCC=gcc
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
CFLAGS= -g -Wall -DVERSION=\"$(GIT_VERSION)\" -std=c++0x #-D_NO_RAZF
INCLUDE=-Isrc/lib -Isrc/bamtools/include -Isrc/bamtools/src -Isrc/ -Isrc/fastahack -Isrc/Complete-Striped-Smith-Waterman-Library/src/ -Isrc/seqan/core/include/ -Isrc/seqan/extras/include
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
FASTA.o:
	cd src/fastahack && make
ssw_cpp.o:
	cd src/Complete-Striped-Smith-Waterman-Library/src && make

SSW = src/Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.o src/Complete-Striped-Smith-Waterman-Library/src/ssw.o
FASTAHACK = src/fastahack/Fasta.o                                                                                                                                                                           
buildWHAMBAM: libbamtools.a FASTA.o ssw_cpp.o
	$(CC) $(CFLAGS) src/lib/*cpp src/bin/multi-wham-testing.cpp $(INCLUDE) $(LIBS) $(FASTAHACK) $(SSW)  -o $(OUTFOLD)WHAM-BAM $(RUNTIME)
buildWHAMBAMD: libbamtools.a FASTA.o ssw_cpp.o
	$(CC) $(CFLAGS) -g -DDEBUG src/lib/*cpp src/bin/multi-wham-testing.cpp $(INCLUDE) $(LIBS) $(FASTAHACK) $(SSW) -o $(OUTFOLD)WHAM-BAM $(RUNTIME)
buildWHAMDUMPER:
	$(CC) $(CFLAGS) -g src/lib/*cpp   src/bin/multi-wham.cpp $(INCLUDE) $(LIBS) -o $(OUTFOLD)WHAM-BAM-DUMPER $(RUNTIME)
buildWHAMBAMGENE:
	$(CC) $(CFLAGS) -g src/lib/*cpp  src/bin/multi-wham-testing-gene.cpp  $(INCLUDE) $(LIBS) -o $(OUTFOLD)WHAM-BAM-GENE $(RUNTIME)
clean:
	-@rm *.a