## Creation of reference genome and alignment of WGS data

### Software utilized:

 - bwa: version .7.7-r441
 - samtools: version 0.1.19-44428cd
 - picard: version 1.138

### Create reference genome

Download the mm9 reference from [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz). Merge it's contents, excluding files containing "random", into a single fasta file, then append the transposon information to that (containined in "transposon.fa"), and finally index with bwa. Below is an exmaple on how to do this:  

```
refdir="/srv/mm9wgs"

#Downloading chromosome reference sequences
curl http://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz > $refdir/chromFa.tar.gz

#Unpacking and merging into standard fasta format
tar -xf $refdir/chromFa.tar.gz -C $refdir
rm $refdir/chr*_random.fa

#Appending transposon sequence
cat $refdur/chr*.fa > $refdir/onc2.fa
cat transposon.fa >> $refdir/onc2.fa

#Indexing reference sequences
bwa index $refdir/onc2.fa
```

### Align data and mark duplicates

Once the reference genome has been created, bwa is used to map reads to the mouse+transposon reference, samtools to sort the mapped 
reads, and picard to dedup sorted reads. Note that you should change the "picarddir" to the directory on your system containing 
picard.jar. Note that you will need to have your paired-end WGS data in a "fastq" directory in your current working directory. 
In the example below, the fastq files are named XXXXX\_R1.fastq and XXXXX\_R2.fastq, where XXXXX is denoted by the "dataset" variable:

```
dataset=Mann_D5SP
picarddir=/nnlab/seqkit/bin
refdir="/srv/mm9wgs"

mkdir -p bam

# Map reads (3 days on 40 processor, 1TB memory rig)
bwa mem -t 40 -M -R "@RG\tID:foo\tSM:bar" $refdir/onc2.fa fastq/$dataset_R1.fastq fastq/$dataset_R2.fastq | samtools view -F4 -bS - > bam/$dataset.mm9.unsorted.bam
samtools index bam/$dataset.mm9.unsorted.bam

# Sort mapped reads (6 hours on 40 processor, 1TB memory rig)
samtools sort -@10 -m 30G -o bam/$dataset.mm9.unsorted.bam - > bam/$dataset.mm9.bam
samtools index bam/$dataset.mm9.bam

# Mark duplicates, a.k.a. dedup (1 day on 40 processor, 1TB memory rig)
java -Djava.awt.headless=true  -Xmx4000m  -jar $picarddir/picard.jar MarkDuplicates I=bam/$dataset.mm9.bam O=bam/$dataset.mm9.dedup.bam M=/dev/null TMP_DIR=/dev/shm
samtools index bam/$dataset.mm9.dedup.bam

```

Due to the large size of each file, it is recommended that you have about 3TB scratch space and 500GB of memory to run the above commands. Note that in the deduping command, we set the tmp directory to "/dev/shm" but you may want to simply use the default setting ("/tmp") depending on your setup or preferences.

