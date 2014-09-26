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
