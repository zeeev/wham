### WHAM-GRAPHENING
======

WHAM-GRAPHENING is under development, but produces more accurate deletion, duplication and inversion calls.

The genotyping module can be slow in high-copy regions.  I recommend filtering SVs prior to genotyping (-k skips genotyping).

Full documentation is coming soon. 

Workflow:

![alt tag](https://github.com/jewmanchue/wham/blob/master/docs/wg.png)

### All WHAM documents can be found on the wiki:

http://jewmanchue.github.io/wham/


### WHAM
======

WHole-genome Alignment Metrics (WHAM) is a structural variant (SV) caller that integrates several sources of mapping information to identify SVs.  WHAM classifies SVs using a flexible and extendable machine-learning algorithm (random forest).  WHAM is not only accurate at identifying SVs, but its association test can identify shared SVs enriched in a cohort of diseased individuals compared to a background of healthy individuals.   

WHAM can be easily run as a stand alone tool or as part of gkno (http://gkno.me) or bcbio-nextgen (http://bcbio-nextgen.readthedocs.org) pipelines.  

Currently we are preparing the WHAM manuscript.  If you would like to use and cite this tool, for the time being, please cite the WHAM github page.


