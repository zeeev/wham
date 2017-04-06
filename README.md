### wham

The wham suite consists of two programs, wham and whamg. wham, the original tool, is a very sensitive method with a high false discovery rate. The second program, whamg, is more accurate and better suited for general structural variant (SV) discovery.
**For most studies, we strongly recommend using whamg**.

Below, we outline the basics of running whamg. **Important sections are highlighted in bold text.** Please cite the [wham paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004572) if you use wham or whamg.

 [![Build Status](https://travis-ci.org/zeeev/wham.svg?branch=master)](https://travis-ci.org/zeeev/wham)  [![Analytics](https://ga-beacon.appspot.com/UA-50580904-5/wham/blob/master/README.md)](https://github.com/zeeev/wham)  [![Code Climate](https://codeclimate.com/github/zeeev/wham/badges/gpa.svg)](https://codeclimate.com/github/zeeev/wham)

### whamg workflow

![alt tag](https://github.com/zeeev/wham/blob/master/docs/github-figure.png)

### Installing whamg

whamg depends on ```CMake``` and ```OpenMP```.  These dependances are commonly preinstalled on various Linux distributions. The following command should install wham into your current working directory:
```
git clone --recursive  https://github.com/zeeev/wham.git; cd wham; make 
```

### Running whamg

whamg uses paired alignments generated with BWA-MEM.  whamg uses the same BAMs as SNV and INDEL calling tools.  Duplicates should be marked or removed, and indel realignment is helpful.  whamg is agnostic regarding the BWA-MEM –M flag (if you don’t know what that means, don’t worry). **It is important that the –R flag in BWA-MEM is used.  whamg requires read-group information. Currently, whamg assumes one library per bam file.**


#### Example 1

```
export EXCLUDE=GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605
whamg -e $EXCLUDE -a Homo_sapiens_assembly19.fasta -f CHM1_1.bam | perl utils/filtWhamG.pl > chm1.vcf  2> chm1.err
```

In the example above, whamg will process all chomosomes in the CHM1 genome. The export statement allows you to mask certain chromosomes/contigs; it may be useful to mask certain chromosomes that are either problematic for read mapping (i.e. rDNA) or short contigs from de novo genome assemblies. The provided example masks human contigs that we have found to be problematic in SV calling. The standard output contains the SV calls and is written to chm1.vcf. The standard error is written to chm1.err, and includes progress updates that can be useful to track runtimes and errors. If you have a trio or a quad and want to joint call, you can pass the bam files as a comma-separated list or a simple text file with the full path to each bam file on each line. **PLEASE NOTE: It is very important that the –e or -c flags are used.** If –e or –c are not specified, there is a chance that whamg will incorrectly estimate the insert sizes of the library.

**A filtering script is in the pipe to increase specificity.**

#### Example 2
```
whamg -x 2 -e $EXCLUDE -a Homo_sapiens_assembly19.fasta -f CHM1_1.bam > chm1.vcf  2> chm1.err
```

In the example above, we are telling whamg to only use two CPUs (-x).

#### Example 3

```
whamg  -g graph.txt -x 2 -e $EXCLUDE -a Homo_sapiens_assembly19.fasta -f CHM1_1.bam > chm1.vcf  2> chm1.err
```

In the example above, each putative structural variant will be written to a flat text file. Individual graphs can be visualized with dotviz or gephi. This output can be helpful for interrogating missing calls or complex structural variants. Do not use this option on a genome-wide run as the file sizes can be extremely large.

### Explanation of options 
 **-f** : A comma-separated list of indexed bam files or a sample text file that looks like:
```
NA12878.sort.bam
NA12776.sort.bam
```
**-a**: The matched reference genome. It is important that the bams were aligned to the same reference sequence.
**-s**: whamg will exit after collecting the basic statistics:
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
**-g**: The file to write the graphs to. Not recommended for whole-genome runs.

**-e**: A comma-separated list of seqids to skip.  

**-c**: A comma-sepearted list of seqids to use for estimating the insert size distribution.

**-r**: Region in ``seqid:start-end`` format. Runs whamg on a specific genomic region.

**-x**: The number of CPUs whamg will attempt to use.  During the first step, whamg reads the entire bam file. If you notice CPU usage dropping, I/O might be swamping out. After the bam files are read, there are several single CPU steps before whamg finishes. The optimal number of CPUs to use really depends on I/O speeds.

**-i**:  whamg uses the BWA-MEM SA tag (default). Older versions of BWA-MEM used a different tag, XP (-i XP).

**-z**:  Sometimes whamg can fail to sample enough reads (low-fold coverage, exome, …). The –z flag forces whamg to keep sampling random regions until it succeeds



### VCF 4.2 output

In this section, each INFO and FORMAT field will be covered. Here is an example whamg VCF header:

```
##fileformat=VCFv4.2
##source=WHAM-GRAPHENING:v1.7.0-225-g1e35-dirty
##reference=Homo_sapiens_assembly19.fasta
##INFO=<ID=A,Number=1,Type=Integer,Description="Total pieces of evidence">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=CF,Number=1,Type=Float,Description="Fraction of reads in graph that cluster with SVTYPE pattern">
##INFO=<ID=CW,Number=5,Type=Float,Description="SVTYPE weight 0-1; DEL,DUP,INV,INS,BND">
##INFO=<ID=D,Number=1,Type=Integer,Description="Number of reads supporting a deletion">
##INFO=<ID=DI,Number=1,Type=Float,Description="Average distance of mates to breakpoint">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=EV,Number=1,Type=Integer,Description="Number everted mate-pairs">
##INFO=<ID=I,Number=1,Type=Integer,Description="Number of reads supporting an insertion">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Number of split-reads supporing SV">
##INFO=<ID=SS,Number=1,Type=Integer,Description="Number of split-reads supporing SV">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=T,Number=1,Type=Integer,Description="Number of reads supporting a BND">
##INFO=<ID=TAGS,Number=.,Type=Integer,Description="SM tags with breakpoint support">
##INFO=<ID=TF,Number=1,Type=Integer,Description="Number of reads mapped too far">
##INFO=<ID=U,Number=1,Type=Integer,Description="Number of reads supporting a duplication">
##INFO=<ID=V,Number=1,Type=Integer,Description="Number of reads supporting an inversion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Per sample SV support">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CHM1
```

| FIELD | DESCRIPTION |
|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| A | Each pair of reads counts as a piece of evidence. The total evidence is summed across the whamg graph structure. Since whamg does not initially cluster similar graph structures, you may see several SVs with low total support that are off by a couple base pairs. Merging SVs is a good idea. |
| CIEND and CIPOS | This field is fixed at [-10, 10]. Most whamg calls are accurate to within a couple base pairs. However, this field is updated during SV merging. |
| CF | Larger SVs (greater than average insert size of the library) generate mate pair patterns. For example, the mates of reads overlapping the 5’ end of a deletion should map downstream of the 3’ breakpoint. Each SV class generates a different pattern. CF is the fraction of the reads in the graph structure that cluster according to SVTYPE. |
| CW | Each mate pair and split read can contribute evidence to one or more SVTYPES. The CW field describes what fraction of the total support belongs to each SVTYPE. In the case of inverted deletions, you might find that both the inversion and deletion have high weights. This field sums to one. |
| D | The number of mate pairs in the graph that support a deletion. |
| DI | This field is related to the CF field. This measures the average distance of the mate pairs to the alternative breakpoint. Inversions can have large DI scores, but deletions should not |
| END | This is a standard VCF field. Insertions start and end at the same position. |
| EV | The number of everted mate pairs (← →). |
| I | The number of reads supporting and inversion. |
| SR | The number of split reads in the graph. This is not necessarily the same as the number of split reads between the breakpoints. |
| SS | The number of reads on the same strand (→ → or ← ←). |
| SVLEN | This is a standard VCF field. |
| SVTYPE | This is a standard VCF field. |
| T | The number of reads supporting and interchomosmal event. |
| TAGS | The individuals (SMs) that are in the graph. |
| TF | The number of mate pairs that map too far away. |
| U | The number of reads that support duplication. |
| V | The number of reads that support an inversion. |
| GT | Genotype. Currently, there is no genotype provided. We are working on a fast and accurate genotyper. |
| DP | Depth. Currently, no depth is provided. |
| SP | This is the number of reads in each individual that supports the exact breakpoint. Because of breakpoint variability, this number might be lower than expected. **Be cautious when filtering on SP.** |


### Filtering

To balance sensitivity and specificity it is wise to filter the raw whamg output.  Filtering always depends on your dataset;  i.e. high depth datasets should be treated differently than low depth datasets.  Similarly, joint calling/family calling should be filtered differently than individual level calling.  Here are some tricks to help guide you.  There are many different tools and strategies that will work for filtering, below are just a few ways you could filter whamg calls.

 

The discussion below covers filtering a high coverage trio, but the concepts will be the same for a single genome.

 

#### Size filtering

 

Whamg emits calls with a wide range of sizes; where larger and smaller events should be subject to more scrutiny.  As a rough rule, we remove calls smaller than 50bp and larger than 2Mbs.  This will discard real calls, but will improve the overall accuracy.

 

#### Total support (INFO field “A”)

 

For a high coverage human genome (~ 50X), filtering calls with less than 5 supporting reads is wise.  Below is a plot of the total support for a joint called trio broken into high and low support bins.  Using a cutoff of 5 removes 20% of the calls (22,955/38,735 remain), which is reasonable given a joint called high coverage trio.

![alt tag](https://github.com/zeeev/wham/blob/master/docs/whamg_total_support_png.png)
 

#### Cross chromosomal mappings

 

Whamg does not call BNDs/translocations, but it is aware of this type of evidence.  We filter out events that have high cross-chromosome mappings.  Pull out the “CW” info field shown below.  In it you will find, listed in order, the fraction of reads supporting  each SV Type.  Remove SVs that have a “BND CW” greater than 0.2.

 

```

##INFO=<ID=CW,Number=5,Type=Float,Description="SVTYPE weight 0-1; DEL,DUP,INV,INS,BND">

```

 

#### Low support across jointly called SVs

 

For joint called variants, low support coming from each individual can sum up to moderate total support.  [We want to remove calls where no individual has significant support?] From the format field, shown below, pull out the “SP” field.  This enumerates the evidence in each individual.  We remove calls where none of the individuals reach an arbitrary level.  This filter is similar to the “Total support filter”

 

```

##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Per sample SV support">
```



---

### wham

wham is the other SV called provided in the suite and represents the codebase published in our PLoS Computational Biology paper. All wham documents can be found on the wiki:

http://zeeev.github.io/wham/

WHole-genome Alignment Metrics (wham) is a structural variant (SV) caller that integrates several sources of mapping information to identify SVs. wham classifies SVs using a flexible and extendable machine-learning algorithm (random forest). wham is not only accurate at identifying SVs, but its association test can identify shared SVs enriched in a cohort of diseased individuals compared to a background of healthy individuals.   

wham can be easily run as a stand-alone tool or as part of gkno (http://gkno.me) or bcbio-nextgen (http://bcbio-nextgen.readthedocs.org) pipelines.
