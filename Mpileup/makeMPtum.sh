#!/bin/sh
samtools mpileup -B -q 1 -f ../Reference/mm9/onc2.fa ../SplitBam/Mann_D5SP/${1}.bam | gzip > ${1}-T.mpileup.gz

