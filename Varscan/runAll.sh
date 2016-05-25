#!/bin/sh
find ../Mpileup/*D5*.mpileup.gz  | perl -i -pe 's/[\/\-]/\n/g' | perl -i -pe 's/mm9/\n/g' | grep dedup | perl -i -pe 's/\.c/c/g' | parallel -j 6 ./runVS.sh

R CMD BATCH dnaCopyAllChroms.R
