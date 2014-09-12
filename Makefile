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
OPENMP=-fopenmp
BOOST=-lboost_program_options
all: createBin buildWHAMBAM buildWHAMBAMGENE

createBin:
	-mkdir bin

buildWHAMBAM:
	$(CC) $(CFLAGS) -g src/lib/*cpp  src/bin/multi-wham-testing.cpp -L $(BAMLIB) -I $(BAMINCLUDE) -I $(WHAMINCLUDE) -lbamtools $(OPENMP) -o $(OUTFOLD)WHAM-BAM $(RUNTIME),src/bamtools/lib
buildWHAMDUMPER:
	$(CC) $(CFLAGS) -g src/lib/*cpp   src/bin/multi-wham.cpp -L $(BAMLIB) -I $(BAMINCLUDE) -I $(WHAMINCLUDE) $(BOOST) -lbamtools $(OPENMP) -o $(OUTFOLD)WHAM-BAM-DUMPER $(RUNTIME),src/bamtools/lib
buildWHAMBAMGENE:
	$(CC) $(CFLAGS) -g src/lib/*cpp  src/bin/multi-wham-testing-gene.cpp -L $(BAMLIB) -I $(BAMINCLUDE) -I $(WHAMINCLUDE) -lbamtools $(OPENMP) -o $(OUTFOLD)WHAM-BAM-GENE $(RUNTIME),src/bamtools/lib
buildWHAMS:
	$(CC) $(CFLAGS) -g src/lib/*cpp  src/bin/single-wham-testing.cpp -L $(BAMLIB) -I $(BAMINCLUDE) -I $(WHAMINCLUDE) -lbamtools $(OPENMP) -o $(OUTFOLD)WHAM-BAM-SINGLE $(RUNTIME),src/bamtools/lib	

runTest:
	bin/WHAM-BAM  -r scaffold974:0-1000000000 -t /home/mshapiro/data/pigeon/bams/large_bams/COLcalRAZDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRARDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAGDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAWDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAQDIAAPE.sort.rmdup.realn.bam  -b /home/mshapiro/data/pigeon/bams/large_bams/COLcalRASDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/PIGtkyRDIDIAAPEI-7.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRADDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAMDIAAPE.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRAODIAAPE_REP.sort.rmdup.realn.bam,/home/mshapiro/data/pigeon/bams/large_bams/COLcalRALDIAAPE.sort.rmdup.realn.bam  > t.txt 2> t.err
#cat scaff79.wham.vcf |  perl -lane '$_ =~/LRT=(.*?);/; print "$F[1]\t$1"' > t2.txt

clean:
	-rm -rf bin/
	-rm t2.txt
	-rm t.txt
	-rm t.err