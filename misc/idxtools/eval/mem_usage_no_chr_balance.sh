#!/bin/bash

#path to the modified minimap2 binary 
MINIMAP=../../minimap2-arm/minimap2

#path to the reference index (download and extract from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz )
REFERENCE=hg38noAlt.fa

#path to the a small read set
READS=../../minimap2-arm/test/q2.fa


###################################################################################

stats_before_opti(){

	TYPE=$1
	PROFILE=map-$TYPE
	FILE_SIZE_LOG=filesize_$TYPE.log
	INDEX_LOG=index_$TYPE.log
	SAM_MAP_LOG=sam_$TYPE.log
	PAF_MAP_LOG=paf_$TYPE.log
	
	echo "filesizes" > $FILE_SIZE_LOG
	echo "index log" > $INDEX_LOG
	echo "Sam map log" > $SAM_MAP_LOG
	echo "Pad map log" > $PAF_MAP_LOG

	#INDEXING
	echo "Generating single reference index"
	/usr/bin/time -v $MINIMAP -a -x $PROFILE $REFERENCE -d hg38noAlt_single.idx > /dev/null 2>> $INDEX_LOG
	ls -l hg38noAlt_single.idx | awk '{print $5}' >> $FILE_SIZE_LOG

	for i in 1700M 850M 350M 140M
	do
		echo "Generating partitioned index for -I $i"
		/usr/bin/time -v $MINIMAP -a -x $PROFILE $REFERENCE -d hg38noAlt_$i.idx -I$i > /dev/null 2>> $INDEX_LOG
		ls -l hg38noAlt_$i.idx | awk '{print $5}' >> $FILE_SIZE_LOG
	done
	#1700M-2parts 850M-4parts 350M-8parts 140M-16parts

	#MAPPING
	echo "Mapping to single REFERENCEerence index : SAM"
	/usr/bin/time -v $MINIMAP -a -x $PROFILE hg38noAlt_single.idx $READS  > groundtruth.sam 2>> $SAM_MAP_LOG

	echo "Mapping to single REFERENCEerence index : PAF"
	/usr/bin/time -v $MINIMAP -x $PROFILE hg38noAlt_single.idx $READS  > groundtruth.paf 2>>  $PAF_MAP_LOG

	for i in 1700M 850M 350M 140M
	do
		echo "Mapping to partitioned index for -I $i with merging: SAM"
		/usr/bin/time -v  $MINIMAP -a -x $PROFILE hg38noAlt_$i.idx $READS  --multi-prefix tmp > merge.$i.sam 2>> $SAM_MAP_LOG
		echo "Mapping to partitioned index for -I $i with merging: PAF"
		/usr/bin/time -v $MINIMAP -x $PROFILE hg38noAlt_$i.idx $READS  --multi-prefix tmp > merge.$i.paf 2>> $PAF_MAP_LOG
	done

	rm *.sam *.paf *.idx
}

logcat(){

	cat paf_pb.log | grep "Maximum resident" | awk '{print $NF/(1024*1024)}' | tr '\n' '\t'
	echo -e -n "\nindex_ont\t"
	cat index_ont.log | grep "Maximum resident" | awk '{print $NF/(1024*1024)}' | tr '\n' '\t'
	echo -e -n "\nsam_ont\t"
	cat sam_ont.log | grep "Maximum resident" | awk '{print $NF/(1024*1024)}' | tr '\n' '\t'
	echo -e -n "\npaf_ont\t"
	cat paf_ont.log | grep "Maximum resident" | awk '{print $NF/(1024*1024)}' | tr '\n' '\t'
	echo ""
}

stats_before_opti "ont"
logcat

