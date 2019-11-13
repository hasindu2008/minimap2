#!/bin/bash

#path to the modified minimap2 binary 
MINIMAP=../minimap2-arm/minimap2

#path to paftools
PAFTOOLS=../minimap2-arm/misc/paftools.js

#path to the k8 javascript shell
K8=../minimap2-arm/misc/k8

#path to the reference index (download and extract from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz )
REFERENCE=hg38noAlt.fa

#path to the synthetic read set (extract pb-1.fa.gz using gunzip)
READS=pb-1.fa


###################################################################################

SINGLE_IDX_AND_PART_IDX_NO_MERGE(){

#generate indexes
echo "Generating single reference index"
$MINIMAP -a -x map-pb $REFERENCE -d hg38noAltpb_single.idx > /dev/null 2> hg38noAltpb_single.idx.log
#1700M-2parts 850M-4parts 350M-8parts 140M-16parts
for i in 1700M 850M 350M 140M
do
	echo "Generating partitioned index for -I $i"
    $MINIMAP -a -x map-pb $REFERENCE -d hg38noAltpb_$i.idx -I$i > /dev/null 2> hg38noAltpb_$i.idx.log
done


#Evaluation of a single reference index
echo "Evaluating SAM output for single reference index"
$MINIMAP -a -x map-pb hg38noAltpb_single.idx $READS -t 8 > groundtruth.sam 2> groundtruth.sam.log
$K8 $PAFTOOLS mapeval groundtruth.sam > groundtruth.sam.eval
echo "Evaluating PAF output for single reference index"
$MINIMAP -x map-pb hg38noAltpb_single.idx $READS -t 8 > groundtruth.paf 2> groundtruth.paf.log
$K8 $PAFTOOLS mapeval groundtruth.paf > groundtruth.paf.eval

#Evaluation of a partitioned index without merging
for i in 1700M 850M 350M 140M
do
	echo "Evaluating partitioned index with no merging for -I $i : SAM"
    $MINIMAP -a -x map-pb hg38noAltpb_$i.idx $READS -t 8 > before.$i.sam 2> before.$i.sam.log
    grep -v "^@" before.$i.sam | sort -k1,1  > sorted.before.$i.sam
    $K8 $PAFTOOLS mapeval sorted.before.$i.sam >  before.$i.sam.eval
	echo "Evaluating partitioned index with no merging for -I $i : PAF"
    $MINIMAP -x map-pb hg38noAltpb_$i.idx $READS -t 8 > before.$i.paf 2> before.$i.paf.log
    sort -k1,1 before.$i.paf > sorted.before.$i.paf
    $K8 $PAFTOOLS mapeval sorted.before.$i.paf >  before.$i.paf.eval
done
}



PART_IDX_MERGED_EVAL(){

#Evaluation of a partitioned index with merging
for i in 1700M 850M 350M 140M
do
	echo "Evaluating partitioned index with merging for -I $i : SAM"
    $MINIMAP -a -x map-pb hg38noAltpb_$i.idx pb-1.fa -t 8 --multi-prefix tmp > merge.$i.sam 2> merge.$i.sam.log
    $K8 $PAFTOOLS mapeval merge.$i.sam >  merge.$i.sam.eval
	echo "Evaluating partitioned index with merging for -I $i : PAF"
    $MINIMAP -x map-pb hg38noAltpb_$i.idx pb-1.fa -t 8 --multi-prefix tmp > merge.$i.paf 2> merge.$i.paf.log
    $K8 $PAFTOOLS mapeval merge.$i.paf >  merge.$i.paf.eval
done
}


SINGLE_IDX_AND_PART_IDX_NO_MERGE
echo "Single reference index and partitioned index without merge has been evaluated"
echo ""
PART_IDX_MERGED_EVAL
echo "Partitioned index with merge has been evaluated"
