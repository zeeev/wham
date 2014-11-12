#### Version 1.0.0

Released October 20th 2014

#####Features:
1. Identification of SV breakpoints.
2. Identification of soft-clipped consensus.
3. Provides the end of breakpoints using supplement alignments.
4. Provides genotypes and genotype likelihoods for multiple individuals.
5. Multi-threading enabled, supported by OpenMP.

#####Known Issues:
1.  For large population calling WHAM can hit the Ulimit.
2.  For deep regions WHAM has considerable slowdown.
3.  Memory becomes an issue for deep regions.

#####Future efforts:
1. Improvements to speed and memory usage.
2. Improve WHAM's break-point joining.
3. Improve the LRT for association testing.


#### Version 1.1.0

Released October 24th 2014

#####Features:
1. WHAM now only reports SVs where the consensus is longer than 10bp.
2. WHAM can now be built in debug mode, which prints all critical data structures
```
make debug
```

#####Known Issues:
Same as Version 1.0.0

#####Future efforts:
Same as Version 1.0.0

#####Bug fixes:
1. Now using discordant mate-pairs correctly.  
2. Major improvements to accuracy on simulated data
3. More SVs have paired breakpoints


#### Version 1.2.0

Released October 31th 2014 - happy Halloween!

#####Features:
1. Increased overall accuracy of SV calls
2. Added attributes to VCF file for SV type classification
3. Added utils/classifier_parse.py for SV classification 
4. Added functionality to readPileUp Class

#####Known Issues:
Same as Version 1.1.0

#####Future efforts:
Same as Version 1.1.0

#####Bug fixes:

#### Version 1.3.0 

Released Nov 4th - happy voting!

#####Features:
1. Increased overall accuracy of SV calls.
   Large calls and translocations require more support.
2. Added another attribute for the classifier.  

#####Future efforts:
1. Validation on published data sets.

#####Bug fixes:

#### Version 1.4.0

Released Nov 12th - happy Veterans Day -1

### Features:
1. Added deltaAfWham for pooled sequencing experiments.
2. Classifier now prints cross-validation information.

#####Future efforts:

1. Adding subsampling to the genotype model.

#####Bug fixes:

1. High depths were causing the genotype likelihood model
to underflow.  This is fixed.  Only the first 1k reads 
are considered for a genotype call.  
2. Fixes to the VCF header again