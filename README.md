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
##INFO=<ID=WAF,Number=3,Type=Float,Description="Allele frequency of: background,target,combined">
##INFO=<ID=GC,Number=2,Type=Integer,Description="Number of called genotypes in: background,target">
##INFO=<ID=AT,Number=15,Type=Float,Description="Attributes for classification">
##INFO=<ID=PU,Number=1,Type=Integer,Description="Number of reads read supporting position">
##INFO=<ID=SU,Number=1,Type=Integer,Description="Number of supplement read supporting position">
##INFO=<ID=CU,Number=1,Type=Integer,Description="Number of neighboring all soft clip clusters across all individuals at pileup position">
##INFO=<ID=SI,Number=1,Type=Float,Description="Shannon entropy">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Number of reads at pileup position across individuals passing filters">
##INFO=<ID=SP,Number=1,Type=String,Description="Support for endpoint;  none:., mp:mate pair, sr:split read">
##INFO=<ID=BE,Number=3,Type=String,Description="Best end position: chr,position,count or none:.">
##INFO=<ID=DI,Number=1,Type=Character,Description="Consensus is from front or back of pileup : f,b">
##INFO=<ID=NC,Number=1,Type=String,Description="Number of soft clipped sequences collapsed into consensus">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT allele">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GL,Number=A,Type=Float,Description="Genotype likelihood">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Number of reads that do not support a SV">
##FORMAT=<ID=NA,Number=1,Type=Integer,Description="Number of reads supporting a SV">
##FORMAT=<ID=NS,Number=1,Type=Integer,Description="Number of reads with a softclip at POS for individual">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Number of reads passing filters">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878.sorted.bam      NA12877.sorted.bam      NA12882.sorted.bam
```

Here is an unclassified wham example:
```
6	
131440530
.	
N	
AGAGTAGAAAGTATCCTTATTTGGTTGGAAAGATCAGATGAAGGGGA	
.
.
LRT=0.3669;WAF=0.250001,0.500001,0.333334;GC=2,1;AT=0.853846,0,0,0,0.146154,0,0.0615385,0.146154,0,0.0846154,0,0,0,0,0,0,0.00769231;SI=1.81279;PU=22;SU=19;CU=140;RD=130;NC=18;SP=sr;BE=6,131440608,19;DI=b;END=131440608;SVLEN=78	
GT:GL:NR:NA:NS:RD	
0/1:-151.971,-24.2602,-331.572:24:11:11:35	0/0:-4.5e-05,-255,-255:45:0:0:45	0/1:-151.971,-29.1122,-428.281:31:11:11:42
```
The first several columns should be familiar to most users as chrom, pos, ID, ref, alt, qual, filter.  The details for the info, format, and genotype fields can be found below.  The genotypes for the trio are 1/0, 0/0 and 0/1. 


====

The info field is comprised of ten key value pairs.

#### LRT:

  This is the likelihood ratio test statistic quantifying the difference in allele frequencies between the cases and controls.
  If you are only using WHAM for its structural breakpoint detection this field will contain a -nan. Both  cases (target) and controls (background) bams must be specified using the -t & -b flags to get an LRT score.
  
#### WAF:

  A comma separated list of allele frequencies of the background, target, combined. If control (background) bams were not specified the allele frequency for the background will be listed as -nan.  WHAM treats each breakpoint as a bi-allelic variant and estimates the frequency based on the genotype counts.
  
#### GC:

  Wham does not call a genotype unless there are 3 reliable reads covering a position.  Genotypes of './.' are no-calls.  CG reports the number of genotypes that could be called reliably.
  
#### AT:

 This field used by wham to determine a classification.  There are 17 read depth normalized counts.  The average user should not worry about this field.
 
#### SI:

  Consensus sequence complexity measured by Shannon's entropy.  Wham removes calls where the entropy is below 0.5.
  
  
#### PU:
  
  The number of primary mappings that support the exact breakpoint reported in the POS field. 
  
#### SU:

  The number of supplementary / split read support for the exact breakpoint reported in the POS field. 
  
#### CU:

  WHAM skips between positions that have soft-clipping.  There are some number of reads that cover a given soft-clipping position.  These reads can have soft-clipping at other locations.  The number of other positions where there are soft-clipping is reported at CU. 
  
#### RD:

  The number of reads covering all individuals joint called
  
  
#### NC:
  
  The number of soft-clipped segments there were collapsed into the consensus sequence.
  
  
#### SP:

  Wham tries to identify the other breakpoint of an SV by first operating on split reads then mate pair mappings.
  If the endpoint was identified by split reads the support tag will be 'sr', if the breakpoint was identified using mate pair mapping the support tag will be 'mp', otherwise the tag will be '.'.  
  
#### BE:

  The best candidate endpoint based on split read clustering.
  
#### DI:
  
  This field denotes is the breakpoint is supported on the 5' of the pielup or 3' end of pileup.  The value 'f' denotes the '5 and 'b' denotes the 3'. 
  
#### END:
  
  The end position of the structural variant.  If there is an interchromosomal SV this field will contain a '.'.
  
#### SVLEN:

  The distance between the two breakpoints wham identified.   If there is an interchromosomal SV or the endpoint could not be identified this field will contain a '.'.
  
#### WC:

  The class of structural variant classified by the random forest.
  
#### WP: 

  The probabilities for each class label generated by the random forest classifier.  
  
====  

The format field is comprised of six colon-delimited fields.

####GT:
  A genotype call.  The nature of the variant is unknown.  WHAM determines the zygosity.
####GL:
  log10 genotype likelihood under a bi-allelic model.
####NR:
  The number of reads that do not support any type of structural variant.
####NA:
  The number of reads that do support any type of structural variant.
####NS:
  The number of primary reads supporting with a soft clip at POS
####RD:
  The number of reads covering the individual at POS after filtering
  




[![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/jewmanchue/wham/trend.png)](https://bitdeli.com/free "Bitdeli Badge")

