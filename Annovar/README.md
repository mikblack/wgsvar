## Varscan directory 

### Files and folders in this directory

 - README.md: this file
 - runAnno.sh: run annovar for a single chromosome to annotate the varscan output (SNPs, InDels and CNV).
 - runAnnoAll.sh: invoke the runAnno.sh script for all chromosomes
 - joinAnnotation.R: R script to combine per-chromosome annotation data, and perform multiple testing adjustment

### Implementation

Only the runAnnoAll.sh scripts needs to be run:

```
nohup ./runAnnoAll.sh &
```

Note that the `joinAnnotation.R` script is specific to the analysis of data from the Sleeping Beauty mouse
leukemia model.  Use in other settings will require modification of the code to match the characteristics 
of the data being analysed. 
