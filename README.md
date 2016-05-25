## WGS variant detection (wgsvar) README

Last modified 2016-05-26

This repository contains a collection of scripts that can be used to perform variant detection in whole genome sequence
data from a pair of tumour and normal samples obtained from a Sleeping Beauty mouse model of cancer.

Github repository: https://github.com/mikblack/wgsvar

### Compute details

_Minimum hardware specs:_

64-bit architecture, 8 GB RAM, 20 GB disk space, Intel Core i5 processor

_Tested software configurations:_

RHEL 6.5, 6.6, 6.7 with gcc/g++ 4.4+ and bash 4.1+

### Dirctory structure

Each subdirectory contains scripts for one aspect of the WGS data analysis workflow.

To step through the workflow, cd into each subdirectory, and run the relevant scripts.  This 
must be done in the correct order, as each step requires files produced by the previous step.

The order is:

 - SplitBam: split BAM files from tumour and normal sample into per-chromosome data
 - Mpileup: generate mpileup files for each chromosome (tumour and normal)
 - Varscan: identify SNVs, InDels and CNV for each chromosome 
 - Annovar: annotate variants identified by Varscan

The content of each directory is described in more detail below.

In addition to the above directories, two additional directories are required:
 
 - The "AnnovarSetup" directory should contain a working installation of the "annovar" software (see below).

 - The "Reference" directory should contain the a copy of the reference genome (in FASTA format) to which the sequence 
   data has been aligned.  For the Sleeping Beauty mouse data, this is the mm9 reference, with the Sleeping Beauty transposon 
   sequence included as an additional chromosome (see file: transposon.fa).
 
This workflow begins with aligned data in BAM format.  An additional directory ("Align") is included here, which contains instructions
for constructing an appropriate reference genome, and performing an alignment of the fastq data.

### Software required

The following software needs to be installed in order to run the analysis workflow.  Versions used for this work are
indicated next to each application.

 - samtools: version 1.1
 - GNU parallel: 20150622
 - varscan: version 2.3.6
 - annovar: downloaded March 2014
 - R: version 3.0.2
 - R packages: dnacopy (Version 1.36.0)

Additional software used for generating alignments:

 - bwa 
 - samtools 
 - picard 

### SplitBam

Within this directory are two subdirectories:

 - Mann_D5SP: holds the spleen (tumour) BAM files, split by chromosome
 - Mann_D6TL: holds the tail (normal) BAM files, split by chromosome

and two scripts:

 - splitBamNorm.sh
 - splitBamTum.sh

Both scripts need to be run:

```
nohup ./splitBamNorm.sh &
nohup ./splitBamTum.sh &
```

Each script points to a single BAM file for each tissue (tumour or normal) - this needs to 
be changed in each script to reflect the location of this file in the filesystem.

The scripts take each BAM file and split them by chromosome to facilitate parallel processing
in subsequent analysis steps.

### Mpileup

This directory conatins the following scripts

 - makeAll.sh: calls the two scripts below to generate mpileups for all normal and tumour samples
 - makeMPnorm.sh: generate mpileup for single chromosome from normal sample
 - makeMPtum.sh: generate mpileup for single chromosome from tumour sample

Only the makeAll.sh script needs to be run:

```
nohup ./makeAll.sh &
```

### Varscan

This directory contains the following scripts/files:

 - runVS.sh: runs varscan application to generate SNP, InDel and Copy Number calls (tumour versus normal)
 - runAll.sh: calls runVS.sh for all chromosomes
 - dnaCopy.R: processes copy number output data from varscan (per chromosome)
 - dnaCopyAllChroms.R: combines per chromosome copy number information 
 - VarScan.v2.3.6.jar: varscan jar file (called within runVS.sh): this needs to be downloaded

To generate the varscan and R output data, only the runAll.sh script needs to be called:

```
nohup ./runAll.sh &
```

### Annovar

This directory contains the following scripts:

 - runAnno.sh: run annovar for a single chromosome to annotate the varscan output (SNPs, InDels and CNV).
 - runAnnoAll.sh: invoke the runAnno.sh script for all chromosomes
 - joinAnnotation.R: R script to combine per-chromosome annotation data, and perform multiple testing adjustment

Only the runAnnoAll.sh scripts needs to be run:

```
nohup ./runAnnoAll.sh &
```

Note that the `joinAnnotation.R` script is specific to the analysis of data from the Sleeping Beauty mouse
leukemia model.  Use in other settings will require modification of the code to match the characteristics 
of the data being analysed. 
