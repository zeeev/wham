The wham suite consists of two programs, wham and whamg.  Wham, has the highest sensitivity, but a much higher false discovery rate.  For general structural variant discovery use whamg. 


###Installing whamg

Whamg depends on CMAKE and OPENMP.  These dependances are commonly found preinstalled on linux systems. 
```
git clone --recursive  https://github.com/zeeev/wham.git ; cd wham ; checkout devel ; make 
```

###Running whamg

Whamg uses paired alignments generated from bwa-mem.  Whamg uses the same bams as SNV and INDEL calling tools.  Duplicates should be marked or removed, and indel realignment is helpful.  Whamg is agnostic regarding the bwa-mem –M flag (if you don’t know what that means don’t worry). 


####Example 

```
export EXCLUDE=GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605

whamg -e $EXCLUDE -a Homo_sapiens_assembly19.fasta –f CHM1_1.bam > chm1.vcf  2> chm1.err

```

In the example above we are running whamg on the whole CHM1 genome.  The standard error is written to chm1.err, this includes progress updates.  The standard output contains the VCF.  If you have a trio or a quad and want to joint call you can pass the bam files as a comma separated list or a simple text file with the full path to each bam file on each line.

#### Example 2
```
whamg -x 2 -e $EXCLUDE -a Homo_sapiens_assembly19.fasta –f CHM1_1.bam > chm1.vcf  2> chm1.err
```

In the example above we are telling whamg to only use two CPUs.

#### Example 3

```
whamg  -g graph.txt -x 2 -e $EXCLUDE -a Homo_sapiens_assembly19.fasta –f CHM1_1.bam > chm1.vcf  2> chm1.err

```

In the example above each putative structural variant will be written to a flat text file.  Individual graphs can be visualized with dotviz or gephi.  This output can be helpful for interrogating missing calls, or complex structural variants.  


### Explanation on options 
 -f : A comma separated list of indexed bam files or a sample text file that looks like:
	```
NA12878.sort.bam
NA12776.sort.bam
```
-a: The matched reference genome.  It is important that the bams were aligned to the same reference sequence.
-s: whamg will exit after collecting the basic statistics:
	```
INFO: for Sample:CHM1
  STATS:    CHM1: mean depth ..............: 38.5064
  STATS:    CHM1: sd depth ................: 8.93088
  STATS:    CHM1: mean insert length: .....: 292.639
  STATS:    CHM1: median insert length ....: 264
  STATS:    CHM1: sd insert length ........: 129.843
  STATS:    CHM1: lower insert cutoff .....: 32.9538
  STATS:    CHM1: upper insert cutoff .....: 617.245
  STATS:    CHM1: average base quality ....: 36.2782
  STATS:    CHM1: average mapping quality .: 57.9975
  STATS:    CHM1: number of reads used ....: 100405
``` 


### wham [![Build Status](https://travis-ci.org/zeeev/wham.svg?branch=master)](https://travis-ci.org/zeeev/wham)

All wham documents can be found on the wiki:

http://zeeev.github.io/wham/

WHole-genome Alignment Metrics (wham) is a structural variant (SV) caller that integrates several sources of mapping information to identify SVs.  Wham classifies SVs using a flexible and extendable machine-learning algorithm (random forest).  WHAM is not only accurate at identifying SVs, but its association test can identify shared SVs enriched in a cohort of diseased individuals compared to a background of healthy individuals.   

WHAM can be easily run as a stand alone tool or as part of gkno (http://gkno.me) or bcbio-nextgen (http://bcbio-nextgen.readthedocs.org) pipelines.  

