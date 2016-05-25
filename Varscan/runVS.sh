#!/bin/sh
# Added awk code to remove bases with 0 coverage in mpileup file, as this
# was causing varscan to crash
gzip -dc ../Mpileup/Mann_D5SP.mm9.${1}-T.mpileup.gz | awk '{if($4 != 0) print $0}' > ${1}-T.mpileup
gzip -dc ../Mpileup/Mann_D6TL.mm9.${1}-N.mpileup.gz | awk '{if($4 != 0) print $0}' > ${1}-N.mpileup

## Define each analysis:

## Somatic variants
# run varscan to produce "standard" otuput
VS1=`echo java -jar VarScan.v2.3.6.jar somatic ${1}-N.mpileup ${1}-T.mpileup ${1} --min-coverage 10 --min-coverage-normal 10 --min-coverage-tumor 10 --min-var-freq 0.1 --min-freq-for-hom 0.75 --somatic-p-value 0.05`
# run varscan to produce vcf output
VS2=`echo java -jar VarScan.v2.3.6.jar somatic ${1}-N.mpileup ${1}-T.mpileup ${1} --min-coverage 10 --min-coverage-normal 10 --min-coverage-tumor 10 --min-var-freq 0.1 --min-freq-for-hom 0.75 --somatic-p-value 0.05 --output-vcf`

## Copy number analysis
VS3=`echo java -jar VarScan.v2.3.6.jar copynumber ${1}-N.mpileup ${1}-T.mpileup ${1}-cnv --min-segment-size 100 --min-coverage 30 --min-base-qual 20 --min-map-qual 20 --max-segment-size 1000`

## Run these analyses in parallel
echo -e $VS1\\n$VS2\\n$VS3 | parallel ::: 

## Complete the copy number analysis
java -jar VarScan.v2.3.6.jar copyCaller ${1}-cnv.copynumber --output-file ${1}-CNV-called --output-homdel-file ${1}-CNV-called-homdel --min-coverage 30 --min-region-size 1000
echo ${1}-CNV-called > input.txt
R CMD BATCH dnaCopy.R

## Clean up
rm ${1}-T.mpileup
rm ${1}-N.mpileup
