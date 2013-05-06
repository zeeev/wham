May 6-th 2013

Purpose:
Wham is a suite of tools designed to help identify and classify alignment regions
that are of poor quality.  Currently there is one program wham.cpp

Usage:

wham mybam.bam scaffoldX 

Output:

A tab delimited text file.
Columns:

1. Seqid.
2. Position in the pileup. 
3. Number of "bad reads" covering the position in the pileup.
4. Number of reads covering the position in the pileup.
5. Fraction of "bad reads"
6. Probability of observing K bad reads out of N trials (See Details)
7. Global fraction of bad reads.

Installing

Dependencies:

Bamtools 

Boost - functions should be stable, but built on 1.42 

g++ -o wham -I ../bamtools-1.0.2/include/ \
-L $PATH/bamtools-1.0.2/lib/ -I ../zk_lib/ -I $PATH/boost/1.42.0/include/ \
-lbamtools wham.cpp ../whamlib/*cpp

Details:

The Samtools flags are hashed across the entire scaffold.  From this the global fraction
of bad reads is calculated.  Then the scaffold is "pileup."  For each position with 
coverage the CDF of the binomial is calculated where:

K = number of bad reads (column 4)
N = number of reads  (column 5)
P = global fraction of bad reads (column 7)

Bad reads include: Mate not mapped, and same strand.  Further information will be added
soon.  IE: the insert size etc....
