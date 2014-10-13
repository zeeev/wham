######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=g++
CFLAGS=-std=c++0x -Wall
INCLUDE=-Isrc/lib -Isrc/bamtools/include -Isrc/bamtools/src -Isrc/seqan/core/include/ -Isrc/seqan/extras/include
OUTFOLD=bin/
LIBS=-L./ -lbamtools -fopenmp -lz -lm
RUNTIME=-Wl,-rpath=src/bamtools/lib/

all: createBin bamtools libbamtools.a buildWHAMBAM clean

createBin:
	-mkdir bin
bamtools:
	cd src/bamtools && mkdir -p build && cd build && cmake .. && make
libbamtools.a: bamtools
	cp src/bamtools/lib/libbamtools.a .
buildWHAMBAM: libbamtools.a
	$(CC) $(CFLAGS) -g src/lib/*cpp  src/bin/multi-wham-testing.cpp $(INCLUDE) $(LIBS)  -o $(OUTFOLD)WHAM-BAM $(RUNTIME)
buildWHAMDUMPER:
	$(CC) $(CFLAGS) -g src/lib/*cpp   src/bin/multi-wham.cpp $(INCLUDE) $(LIBS) -o $(OUTFOLD)WHAM-BAM-DUMPER $(RUNTIME)
buildWHAMBAMGENE:
	$(CC) $(CFLAGS) -g src/lib/*cpp  src/bin/multi-wham-testing-gene.cpp  $(INCLUDE) $(LIBS) -o $(OUTFOLD)WHAM-BAM-GENE $(RUNTIME)


runTest:
	bin/WHAM-BAM  -r scaffold974:0-1000000000 -t /home/mshapiro/data/pigeon/bams/large_bams/COLcalRAZDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRARDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAGDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAWDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAQDIAAPE.sort.rmdup.realn.bam  -b /home/mshapiro/data/pigeon/bams/large_bams/COLcalRASDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/PIGtkyRDIDIAAPEI-7.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRADDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAMDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAODIAAPE_REP.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRALDIAAPE.sort.rmdup.realn.bam  > t.txt 2> t.err
#cat scaff79.wham.vcf |  perl -lane '$_ =~/LRT=(.*?);/; print "$F[1]\t$1"' > t2.txt

clean:

	-@rm t2.txt
	-@rm t.txt
	-@rm t.err
	-@rm *.a
