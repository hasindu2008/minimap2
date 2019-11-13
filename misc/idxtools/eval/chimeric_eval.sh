#!/bin/bash

fastq=bigRead.fq
name=${fastq%%.*}

#path to the modified minimap2 binary
MINIMAP=../minimap2-arm/minimap2

#path to paftools
PAFTOOLS=../minimap2-arm/misc/paftools.js

#path to the reference index (download and extract from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz )
REFERENCE=hg38noAlt.fa


###########################################################################################

SINGLE_IDX_AND_PART_IDX_NO_MERGE(){

#generate indexes
echo "Generating single reference index"
$MINIMAP -a -x map-ont $REFERENCE -d hg38noAlt_single.idx > /dev/null 2> hg38noAlt_single.idx.log
#1700M-2parts 850M-4parts 350M-8parts 140M-16parts
for i in 1700M 850M 350M 140M
do
	echo "Generating partitioned index for -I $i"
    $MINIMAP -a -x map-ont $REFERENCE -d hg38noAlt_$i.idx -I$i > /dev/null 2> hg38noAlt_$i.idx.log
done


#Evaluation of a partitioned index without merging
echo "Running for single reference index"
$MINIMAP -x map-ont hg38noAlt_single.idx $fastq -t 1 > groundtruth.paf 2> groundtruth.paf.log

#split index results before modifications to minimap
for i in 1700M 850M 350M 140M
do
	echo "Running for partitioned index with no merging for -I $i"
    $MINIMAP -x map-ont hg38noAlt_$i.idx $fastq -t 1 > before.$i.paf 2> before.$i.paf.log
    sort -k1,1 before.$i.paf > sorted.before.$i.paf
done
}


PART_IDX_MERGED_EVAL(){

#Evaluation of a partitioned index with merging
for i in 1700M 850M 350M 140M
do
	echo "Running for partitioned index with merging for -I $i"
    $MINIMAP -x map-ont hg38noAlt_$i.idx $fastq  -t 1 --multi-prefix tmp > merge.$i.paf 2> merge.$i.paf.log
done
}

SINGLE_IDX_AND_PART_IDX_NO_MERGE
echo "Single reference index and partitioned index without merge done"
echo ""
PART_IDX_MERGED_EVAL
echo "Partitioned index with merge done"


