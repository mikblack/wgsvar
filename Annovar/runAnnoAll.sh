#!/bin/sh
find "../Varscan/" -name "*.vcf" | perl -i -pe 's/\.\.\/Varscan\///g' | parallel -j 10 ./runAnno.sh
