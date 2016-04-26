######################################

# Makefile written by Zev Kronenberg #

#     zev.kronenberg@gmail.com       #

######################################



CC=gcc
CXX=g++
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
CFLAGS=  -fstack-protector-all -Wall -DVERSION=\"$(GIT_VERSION)\" -DFAST -std=c++0x  -Wno-sign-compare
INCLUDE=-Isrc/lib -Isrc/bamtools/include -Isrc/bamtools/src -Isrc/ -Isrc/fastahack -Isrc/Complete-Striped-Smith-Waterman-Library/src/ -Isrc/seqan/core/include/ -Isrc/seqan/extras/include
OUTFOLD=bin/
LIBS=-L./  -fopenmp -lz -lm
RUNTIME=-Wl,-rpath=src/bamtools/lib/



all: mvSSW createBin bamtools libbamtools.a buildWHAMBAM whamGraph buildMerge clean
debug: mvSSW createBin bamtools libbamtools.a buildWHAMBAMD graphDebug buildMerge clean

mvSSW:
	cp src/lib/ssw.c src/Complete-Striped-Smith-Waterman-Library/src
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
	$(CXX) $(CFLAGS) src/lib/*cpp src/bin/multi-wham-testing.cpp $(INCLUDE) $(LIBS) $(FASTAHACK) $(SSW)  -o $(OUTFOLD)WHAM-BAM $(RUNTIME)
buildWHAMBAMD: libbamtools.a FASTA.o ssw_cpp.o
	$(CXX) $(CFLAGS) -g -DDEBUG src/lib/*cpp src/bin/multi-wham-testing.cpp $(INCLUDE) $(LIBS) $(FASTAHACK) $(SSW) -o $(OUTFOLD)WHAM-BAM $(RUNTIME)
buildWHAMDUMPER:
	$(CXX) $(CFLAGS) -g src/lib/*cpp   src/bin/multi-wham.cpp $(INCLUDE) $(LIBS) -o $(OUTFOLD)WHAM-BAM-DUMPER $(RUNTIME)
buildWHAMBAMGENE:
	$(CXX) $(CFLAGS) -g src/lib/*cpp  src/bin/multi-wham-testing-gene.cpp  $(INCLUDE) $(LIBS) -o $(OUTFOLD)WHAM-BAM-GENE $(RUNTIME)
whamGraph:
	$(CXX) $(CFLAGS) -O3 src/lib/*cpp src/bin/graph-er.cpp $(INCLUDE) $(LIBS) $(FASTAHACK) $(SSW)  -o $(OUTFOLD)WHAM-GRAPHENING $(RUNTIME)
graphDebug:
	$(CXX) $(CFLAGS) src/lib/*cpp src/bin/graph-er.cpp src/bamtools/lib/libbamtools.a $(INCLUDE) $(LIBS) $(FASTAHACK) $(SSW)  -o $(OUTFOLD)WHAM-GRAPHENING $(RUNTIME) -g

buildTest:
	$(CXX) -g -I/home/zkronenb/tools/gtest-1.7.0/include/ -L/home/zkronenb/tools/gtest-1.7.0/build -
buildMerge:
	$(CXX) $(CFLAGS) $(INCLUDE) $(LIBS) src/bin/mergeIndv.cpp src/lib/split.cpp -o $(OUTFOLD)mergeIndvs


clean:
	-@rm *.a