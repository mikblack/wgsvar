#!/bin/sh

mkdir -p Mann_D6TL
cd Mann_D6TL

bamFile=../Align/bam/Mann_D6TL.mm9.dedup.bam

for i in `seq 1 19`;
do
   samtools view -bh $bamFile chr${i} > Mann_D6TL.mm9.chr${i}.dedup.bam
   samtools index Mann_D6TL.mm9.chr${i}.dedup.bam
done    
samtools view -bh $bamFile chrX > Mann_D6TL.mm9.chrX.dedup.bam
samtools index Mann_D6TL.mm9.chrX.dedup.bam
samtools view -bh $bamFile chrY > Mann_D6TL.mm9.chrY.dedup.bam
samtools index Mann_D6TL.mm9.chrY.dedup.bam

