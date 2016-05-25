## SplitBam directory 

### Files and folders in this directory

 - Mann_D5SP: holds the spleen (tumour) BAM files, split by chromosome
 - Mann_D6TL: holds the tail (normal) BAM files, split by chromosome
 - README.md: this file
 - splitBamNorm.sh: script to split normal samples into per-chromosome BAM files
 - splitBamTum.sh: script to split tumour samples into per-chromosome BAM files

### Implementation

Both scripts need to be run:

```
nohup ./splitBamNorm.sh &
nohup ./splitBamTum.sh &
```

Each script points to a single BAM file for each tissue (tumour or normal) - this needs to 
be changed in each script to reflect the location of this file in the filesystem.

