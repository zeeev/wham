Purpose:
------
WHole-genome Aligment Metricts (wham BAM) is a suite of tools designed to identify anomalies in binary alignment/mapping files (BAM).  

The toolkit has three programs:

### wham 

Wham is designed to look for anomalies within a single bam file

### whamy

Whamy is designed to look for anomalies across groups of bams.  

### raw

the latest tool that incoperates a statistical model

Installing:
-----

###Dependencies:

  Cmake:
    Cmake is designed to make compiling code easy across different operation systems. It can be installed from yum or other repository managers.  To check if cmake is installed type cmake --version.  Cmake 2.6 or higher is required.
    

  Boost:
    Boost is used for it's program opitions and statistical distributions.  If boost is not installed as root. you might need to create an enviromental variable for boost:
    
    bash :
    
    export BOOST_ROOT=/xxx/boost/
    
    tcsh :
    
    setenv BOOST_ROOT /xxx/boost/
    
  LibZ:
    


###Bamtools:

1.  git clone https://github.com/pezmaster31/bamtools
2.  cd bamtools
3.  mkdir build
4.  cd build
5.  cmake ..
6.  sudo make install


For non-root users:

bash users:

add the following to your bash profile

    export PATH=/xxx/bamtools/include:$APTH
    export CPATH=/xxx/bamtools/lib:$CPATH
    export LD_LIBRARY_PATH=/xxx/bamtools/lib:$LD_LIBRARY_PATH


tcsh user (similar to C shell):
 
add the following to your profile:  

    setenv CPATH           {$CPATH}:/xxx/bamtools/include
    setenv  LIBRARY_PATH    {$LIBRARY_PATH}:/xxx/bamtools/lib
    setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:/xxx/bamtools/lib 

###wham BAM:
   
1.  git clone git://github.com/jewmanchue/wham.git
2.  cd wham
3.  mkdir build
4.  cd build
5.  cmake ..
6.  make 
7.  cd ..
8.  rm -rf build

Testing:
-----

sh /test/run-test.sh

Basic usage:
-----

### wham

    wham mybam.bam mybam.bai scaffoldX > stdout 2> stderr

### whamy

Whamy requires sorted BAMs.  To improve whamy's speeds it is a good idea to remove duplicates and index the BAMS (not required).  If the scaffold flag 
is not set whamy will run the whole genome.

    whamy --target a.bam,b.bam,c.bam --background d.bam,e.bam,f.bam --scaffold chr1 > stdout 2> stderr

Output:
-----

###wham

A tab delimited text file.
Columns:

1. Seqid.
2. Position. 
3. Number of "bad reads" covering the position in the pileup.
4. Number of reads covering the position in the pileup.
5. Fraction of "bad reads"
6. Probability of observing K bad reads out of N trials (See Details)
7. Global fraction of bad reads.

###whamy

A tab delimited text file.
Columns:

1. Seqid.
2. Position. 
3. Number of unmapped mates (target)
4. Number of unmapped mates (background)
5. Number of same strand mates (target)
6. Number of same strand mates (background)
7. Number of cross seqid mapped mates (target)
8. Number of cross seqid mapped mates (background)
9. Depth - Number of reads covering position (target)
10. Depth - Number of reads covering position (background)
11. Average insert length (target)
12. Average insert length (background)
13. Average mapping quality (target)
14. Average mapping quality (background)
15. Euclidean distances based on columns (3,4;5,6;7,8;13,14)


Details:
-----

###wham

The Samtools flags are hashed across the entire scaffold.  From this the global fraction
of bad reads is calculated.  Then the scaffold is "pileup."  For each position with 
coverage the CDF of the binomial is calculated where:

K = number of bad reads (column 4)
N = number of reads  (column 5)
P = global fraction of bad reads (column 7)

Bad reads include: Mate not mapped, and same strand.  Further information will be added
soon.  IE: the insert size etc....
