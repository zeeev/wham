######################################
# Makefile written by Zev Kronenberg #
#     zev.kronenberg@gmail.com       #
######################################

CC=g++
CFLAGS=-std=c++0x
BAMINCLUDE=src/bamtools/include/
BAMLIB=src/bamtools/lib/
WHAMINCLUDE=src/lib
OUTFOLD=bin/
RUNTIME= -Wl,-rpath

all: createBin buildWHAMBAM

createBin:
	mkdir bin

buildWHAMBAM:
	c++ $(CFLAGS) -g src/lib/*cpp  src/bin/multi-wham-testing.cpp -L $(BAMLIB) -I $(BAMINCLUDE) -I $(WHAMINCLUDE) -lbamtools -o $(OUTFOLD)WHAM-BAM $(RUNTIME),src/bamtools/lib



#	./a.out  -t /home/mshapiro/data/pigeon/bams/large_bams/COLcalRARDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAGDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAWDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAQDIAAPE.sort.rmdup.realn.bam  -b /home/mshapiro/data/pigeon/bams/large_bams/COLcalRASDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/PIGtkyRDIDIAAPEI-7.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRADDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAMDIAAPE.sort.rmdup.realn.bam  > t.txt 2> t.err
#perl -lane '$_ =~ /RA=(.*?);/; print $1' 