## Varscan directory 

### Files and folders in this directory

 - README.md: this file
 - runVS.sh: runs varscan application to generate SNP, InDel and Copy Number calls (tumour versus normal)
 - runAll.sh: calls runVS.sh for all chromosomes
 - dnaCopy.R: processes copy number output data from varscan (per chromosome)
 - dnaCopyAllChroms.R: combines per chromosome copy number information 
 - VarScan.v2.3.6.jar: varscan jar file (called within runVS.sh): this needs to be downloaded

### Implementation

To generate the varscan and R output data, only the runAll.sh script needs to be called:

```
nohup ./runAll.sh &
```

