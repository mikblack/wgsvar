#!/bin/sh

find ../SplitBam/Mann_D5SP -name "*.bam" | perl -i -pe 's/\.bam//g' | perl -i -pe 's/\//\n/g' | grep dedup | parallel -j 6 ./makeMPtum.sh

find ../SplitBam/Mann_D6TL -name "*.bam" | perl -i -pe 's/\.bam//g' | perl -i -pe 's/\//\n/g' | grep dedup | parallel -j 6 ./makeMPnorm.sh

