#!/bin/sh

mkdir -p Mann_D5SP
cd Mann_D5SP

bamFile=../Align/bam/Mann_D5SP.mm9.dedup.bam

for i in `seq 1 19`;
do
   samtools view -bh $bamFile chr${i} > Mann_D5SP.mm9.chr${i}.dedup.bam
   samtools index Mann_D5SP.mm9.chr${i}.dedup.bam
done    
samtools view -bh $bamFile chrX > Mann_D5SP.mm9.chrX.dedup.bam
samtools index Mann_D5SP.mm9.chrX.dedup.bam
samtools view -bh $bamFile chrY > Mann_D5SP.mm9.chrY.dedup.bam
samtools index Mann_D5SP.mm9.chrY.dedup.bam

