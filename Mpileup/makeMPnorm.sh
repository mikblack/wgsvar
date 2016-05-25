#!/bin/sh
samtools mpileup -B -q 1 -f ../Reference/mm9/onc2.fa ../SplitBam/Mann_D6TL/${1}.bam | gzip > ${1}-N.mpileup.gz

