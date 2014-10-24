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
2. Improve WHAM's breakpoint joining.
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
