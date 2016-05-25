#!/bin/sh
sample=`echo $1 | perl -i -pe 's/.vcf//g'`
echo $sample

grep SOMATIC ../Varscan/${sample}.vcf > ${sample}.SOMATIC.vcf
../AnnovarSetup/convert2annovar.pl -format vcf4old ${sample}.SOMATIC.vcf > ${sample}.SOMATIC.anno
../AnnovarSetup/annotate_variation.pl --buildver mm9 ${sample}.SOMATIC.anno ../AnnovarSetup/mousedb/

