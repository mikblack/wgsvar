## Mpileup directory 

### Files and folders in this directory

 - README.md: this file
 - makeAll.sh: calls the two scripts below to generate mpileups for all normal and tumour samples
 - makeMPnorm.sh: generate mpileup for single chromosome from normal sample
 - makeMPtum.sh: generate mpileup for single chromosome from tumour sample

### Implementation

Only the makeAll.sh script needs to be run:

```
nohup ./makeAll.sh &
```

