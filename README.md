### What is WHAM?

WHole-genome Alignment Metrics (WHAM) is a structural variant (SV) caller that integrates several sources of mapping information to identify SVs.  WHAM classifies SVs using a flexible and extendable machine-learning algorithm (random forest).  WHAM is not only accurate at identifying SVs, but its association test can identify shared SVs enriched in a cohort of diseased individuals compared to a background of healthy individuals.   

WHAM can be easily run as a stand alone tool or as part of gkno (http://gkno.me) or bcbio-nextgen (http://bcbio-nextgen.readthedocs.org) pipelines.  

Currently we are preparing the WHAM manuscript.  If you would like to use and cite this tool, for the time being, please cite the WHAM github page.

### Installing WHAM

WHAM builds in two simple steps... assuming that you have the dependencies.  Most linux environments have cmake, and OpenMP the two requirements for WHAM.  If an install fails you might have to install a few libraries. 

#### Dependencies

For classification of SV variant types, a python distribution with the scikit-learn and argparse packages are rquired. These come standard with the Anaconda release of python, which can be obtained here:

http://continuum.io/downloads

If you want to continue with your current python distribution use the following commands (NOTE: requires Python (>= 2.6 or >= 3.3), NumPy (>= 1.6.1), SciPy (>= 0.9) :

```
pip install -U numpy scipy scikit-learn argparse
```

After dependency installation, clone the git repository with the following commands:

```
git clone --recursive git@github.com:jewmanchue/wham.git
cd wham
git checkout v1.5.0
make
```
If you get an error while downloading wham try https:

```
git clone --recursive https://github.com/jewmanchue/wham.git
cd wham
git checkout v1.5.0
make
```

To make sure you checked out the correct branch you should see v1.0.0 by typing:

```
git name-rev --tags --name-only $(git rev-parse HEAD)
```


### WHAM input

WHAM uses soft clipping information and supplementary alignments provided by mapping software. 
We recommend that you use BWA-mem for accurate structural variant calls.  BAM files must be
sorted, duplicates marked or removed and the BAM files need to be indexed.  The sorted
BAM files need to have the HD:SO tag, this is provided by samtools 0.19 or later.  

### Running WHAM

usage statement:

```
usage  : WHAM-BAM -x <INT> -r <STRING>     -e <STRING>  -t <STRING>    -b <STRING>

example: WHAM-BAM -x 20    -r chr1:0-10000 -e genes.bed -t a.bam,b.bam -b c.bam,d.bam

required: t <STRING> -- comma separated list of target bam files
option  : b <STRING> -- comma separated list of background bam files
option  : r <STRING> -- a genomic region in the format "seqid:start-end"
option  : x <INT>    -- set the number of threads, otherwise max
option  : e <STRING> -- a bedfile that defines regions to score

Version 0.0.1 ; Zev Kronenberg; zev.kronenberg@gmail.com
```

#### Running on test data

We have supplied a test dataset for you to test your installation and get familiar with usage of WHAM. Try the following commands, run from the wham directory. Note that the commands scroll horizontally. 

```
./bin/WHAM-BAM -t test/simulations/inv-5xsimA.sort.rmdup.chr22.bam > test/simulations/inv-5xsimA.sort.rmdup.chr22.test.vcf 
python utils/classify_WHAM_vcf.py test/simulations/inv-5xsimA.sort.rmdup.chr22.test.vcf utils/WHAM_training_data.txt > test/simulations/inv-5xsimA.sort.rmdup.chr22.CALLS.vcf
```
The first command runs WHAM on your input bam file and outputs an unsorted VCF file. The second command runs the WHAM classifier, which uses information from the SV calls and from a training dataset to classify each SV call as a specific SV type (insertions, deletions, etc.) that can be retrieved from the VCF file with a new key 'WC' in the info field.

### Understanding WHAM output

WHAM outputs an unsorted VCFv4.1 file.  Below is an example header.  

```
##fileformat=VCFv4.1
##INFO=<ID=LRT,Number=1,Type=Float,Description="Likelihood Ratio Test Statistic">
##INFO=<ID=AF,Number=3,Type=Float,Description="Allele frequency of: background,target,combined">
##INFO=<ID=GC,Number=2,Type=Integer,Description="Number of called genotypes in: background,target">
##INFO=<ID=CU,Number=1,Type=Integer,Description="Number of neighboring soft clip clusters across all individuals at pileup position ">
##INFO=<ID=ED,Number=.,Type=String,Description="Colon separated list of potenial paired breakpoints, in the format: seqid,pos,count">
##INFO=<ID=BE,Number=2,Type=String,Description="Best end position: chr,position,count">
##INFO=<ID=NC,Number=1,Type=String,Description="Number of soft clipped sequences collapsed into consensus">
##INFO=<ID=AT,Number=13,Type=Float,Description="Attributes for classification">
##INFO=<ID=WC,Number=1,Type=String,Description="WHAM classifier vairant type">
##INFO=<ID=WP,Number=4,Type=Float,Description="WHAM probability estimate for each structural variant classification from RandomForest model">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Pseudo genotype">
##FORMAT=<ID=GL,Number=A,Type=Float,Description="Genotype likelihood ">
##FORMAT=<ID=FR,Number=1,Type=Float,Description="Fraction of reads with soft or hard clipping">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Number of reads supporting a SV">
##FORMAT=<ID=NA,Number=1,Type=Integer,Description="Number of reads that do not support a SV">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of reads with mapping quality greater than 0">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878.sorted.bam      NA12877.sorted.bam NA12882.sorted.bam
```

Here is an example SV identified by WHAM:
```
13
51069351
.
N
CTTAGAAATTTCTTCCACCAGATTCCCTNGAACTTGNCATCNCTCTTAAGT
.
.
LRT=nan;AF=0.250001,1,0.500001;GC=2,1;AT=0.671756,0,0,0,0.274809,0.0152672,0.183206,0.274809,0,0.0916031,0.0534351,0,0,0,0;CU=125;NC=35;ED=13,51075079,36:13,51075111,1:;BE=13,51075079,36;WC=DUP;WP=0.0,0.976,0.02,0.004
GT:GL:NR:NA:DP:FR
1/1:-255,-255,-0.693163:0:16:17:1
0/0:-27.6311,-32.5779,-621.698:45:2:47:0.0425532
0/1:-303.941,-34.6574,-386.834:28:22:50:0.44
```
The first several columns should be familiar to most users as chrom, pos, ID, ref, alt, qual, filter.  The details for the info, format, and genotype fields can be found below. In this example WHAM has annotated a  a duplication (info: WC=DUP) with a probability of 0.976 (WP=0.0,0.976,0.02,0.004) starting at 51,069,351 and ending at 51,075,079 on chromosome 13.  The genotypes for the trio are 1/1, 0/0 and 0/1. 


====

The info field is comprised of ten key value pairs.

#### LRT:

  This is the likelihood ratio test statistic quantifying the difference in allele frequencies between the cases and controls.
  If you are only using WHAM for its structural breakpoint detection this field will contain a -nan. Both  cases (target) and controls (background) bams must be specified using the -t & -b flags to get an LRT score.
  
#### AF:

  A comma separated list of allele frequencies of the background, target, combined. If control (background) bams were not specified the allele frequency for the background will be listed as -nan.  WHAM treats each breakpoint as a bi-allelic variant and estimates the frequency based on the genotype counts.
  
#### GC:

  WHAM does not call a genotype unless there are 3 reliable reads covering a position.  Genotypes of './.' are no-calls.  CG reports the number of genotypes that could be called reliably.
  
#### CU:

  WHAM skips between positions that have soft-clipping.  There are some number of reads that cover a given soft-clipping position.  These reads can have soft-clipping at other locations.  The number of other positions where there are soft-clipping is reported at CU. 

#### NC:
  
  The number of soft-clipped segments there were collapsed into the consensus sequence.
  
#### ED:
  
  All the split read clusters that could be the end position of the SV.  There are three pieces of information in this field: seqid, pos, and the number of split reads supporting this position.  
  
#### BE:

  The best candidate endpoint based on split read clustering.
  
#### WC:

  The class of structural variant classified by the random forest.
  
#### WP: 

  The probabilities for each class label generated by the random forest classifier.  
  
====  

The format field is comprised of 7 colon-delimited fields.

####GT:
  A genotype call.  The nature of the variant is unknown.  WHAM determines the zygosity.
####GL:
  log10 genotype likelihood under a bi-allelic model.
####FR:
  Fraction of the primary alignments that have soft-clipping.
####NR:
  The number of reads that do not support any type of structural variant.
####NA:
  The number of reads that do support any type of structural variant.
####DP:
  The number of high-quality primary reads covering the position.

