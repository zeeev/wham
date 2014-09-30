### WHAM

Structural variants (SVs) have largely been ignored in Genome Wide Association Studies.  SVs are difficult to call from NGS data.  Often the same SV allele is assigned different breakpoints across individuals.  WHAM takes a different approach by directly applying association tests to BAM files.  WHAM does not call SVs, but rather examines the patterns of mate pair mapping across cases and controls.  If there is an enrichment of reads supporting a SV in the cases or controls WHAM performs an association test.


### INSTALLING WHAM

WHAM builds in two simple steps... assuming that you have the depedancies.  Most linux enviroments have cmake, and OpenMP the two requiresments for WHAM.  If an install fails you might have to install a few libraries. 

```
git clone --recursive https://github.com/jewmanchue/wham.git wham
cd wham
make
```

### RUNNING WHAM

usage statement:

```
usage: WHAM -t <STRING> -b <STRING> -r <STRING>

option: t <STRING> -- comma separated list of target bam files
option: b <STRING> -- comma separated list of background bam files
option: r <STRING> -- a genomic region in the format "seqid:start-end"

Version 0.0.1 ; Zev Kronenberg; zev.kronenberg@gmail.com
```

### UNDERSTANDING WHAM OUTPUT

WHAM outputs an unsorted VCFv4.1 file.  Below is an example header.  

```
##fileformat=VCFv4.1
##INFO=<ID=LRT,Number=1,Type=Float,Description="Likelihood Ratio Test Statistic">
##INFO=<ID=AF,Number=3,Type=Float,Description="Allele frequency of: background,target,combined">
##INFO=<ID=GC,Number=2,Type=Integer,Description="Number of called genotypes in: background,target">
##INFO=<ID=CU,Number=1,Type=Integer,Description="Number of neighboring soft clip clusters across all individuals at pileup position ">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Pseudo genotype">
##FORMAT=<ID=GL,Number=A,Type=Float,Desciption="Genotype likelihood ">
##FORMAT=<ID=FR,Number=1,Type=Float,Description="Fraction of reads with soft or hard clipping">
##FORMAT=<ID=NR,Number=1,Type=Integer,Description="Number of reads supporting no SV">
##FORMAT=<ID=NA,Number=1,Type=Integer,Description="Number of reads supporting no SV">
##FORMAT=<ID=CL,Number=1,Type=Integer,Description="Number of bases that have been soft clipped">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of reads with mapping quality greater than 0">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878D_HiSeqX_R1.bam
```

#### LRT:

  This is the likelihood ratio test statistic quantifying the difference in allele frequencies between the cases and controls.
  If you are only using WHAM for its structural breakpoint detection this field will contain a -nan. Both  casess (target) and controls (background) bams must be specified using the -t & -b flags to get an LRT score.
  
#### AF:

  A comma seperated list of allele frequences of the background, target, combined. If control (background) bams were not specified the allele frequency for the background will be listed as -nan.  WHAM treats each breakpoin as a bi-allelic variant and esimates the frequency based on the genotype counts.
  
### GC:

  WHAM does not call a genotype unless there are 3 reliable reads covering a position.  Genotypes of './.' are no-calls.  CG reports the number of genotypes that could be called reliably.
  
### CU

  WHAM skips between positions that have soft-clipping.  There are some number of reads that cover a given soft-clipping position.  These reads can have soft-clipping at other locations.  The number of other positions where there are soft-clipping is reported at CU. 

  
