### WHAM-BAM

Structural variants (SVs) have largely been ignored in Genome Wide Association Studies.  SVs are difficult to call from NGS data.  Often the same SV allele is assigned different breakpoints across individuals.  WHAM-BAM takes a different approach by directly applying association tests to BAM files.  WHAM-BAM does not call SVs, but rather examines the patterns of mate pair mapping across cases and controls.  If there is an enrichment of reads supporting a SV in the cases or controls WHAM-BAM performs an association test.


### INSTALLING WHAM-BAM

WHAM-BAM builds in two simple steps... assuming that you have the depedancies.  Most linux enviroments have cmake, and OpenMP the two requiresments for WHAM-BAM.  If an install fails you might have to install a few libraries. 

```
git clone --recursive git@github.com:jewmanchue/wham.git wham
cd wham
make
```

### RUNNING WHAM-BAM

usage statement:

```
usage: WHAM-BAM -t <STRING> -b <STRING> -r <STRING>

option: t <STRING> -- comma separated list of target bam files
option: b <STRING> -- comma separated list of background bam files
option: r <STRING> -- a genomic region in the format "seqid:start-end"

Version 0.0.1 ; Zev Kronenberg; zev.kronenberg@gmail.com
```
