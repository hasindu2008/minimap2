#!/bin/bash

#path to the modified minimap2 binary 
MINIMAP=../../minimap2-arm/minimap2

#path to chromosome balancing indexer
DIVIDER=../../minimap2-arm/misc/idxtools/divide_and_index.sh

#path to the reference index (download and extract from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz )
REFERENCE=hg38noAlt.fa

#path to the a small read set
READS=../..//minimap2-arm/test/q2.fa


###################################################################################


stats_after_opti(){
	
	TYPE=$1
	PROFILE=map-$TYPE
	FILE_SIZE_LOG=filesize_$TYPE.log
	INDEX_LOG=index_$TYPE.log
	SAM_MAP_LOG=sam_$TYPE.log
	PAF_MAP_LOG=paf_$TYPE.log
	
	
	echo "filesizes" > $FILE_SIZE_LOG
	echo "index log" > $INDEX_LOG
	echo "Sam map log" > $SAM_MAP_LOG
	echo "Paf map log" > $PAF_MAP_LOG

	#INDEXING
	echo "Generating single reference index"
	/usr/bin/time -v $MINIMAP -a -x $PROFILE $REFERENCE -d hg38noAlt_single.idx > /dev/null 2>> $INDEX_LOG
	ls -l hg38noAlt_single.idx | awk '{print $5}' >> $FILE_SIZE_LOG

	for i in 2 4 8 16
	do
		echo "Generating partitioned index for $i parts"
		/usr/bin/time -v $DIVIDER $REFERENCE $i hg38noAlt_$i.idx $MINIMAP $PROFILE  2>> $INDEX_LOG 
		ls -l hg38noAlt_$i.idx | awk '{print $5}' >> $FILE_SIZE_LOG
	done

	#MAPPING
	echo "Mapping to single reference index : SAM"
	/usr/bin/time -v $MINIMAP -a -x $PROFILE hg38noAlt_single.idx $READS  > groundtruth.sam 2>> $SAM_MAP_LOG

	echo "Mapping to single reference index : PAF"
	/usr/bin/time -v $MINIMAP -x $PROFILE hg38noAlt_single.idx $READS  > groundtruth.paf 2>>  $PAF_MAP_LOG

	for i in 2 4 8 16
	do
		echo "Mapping to partitioned index  for $i parts with merging: SAM"
		/usr/bin/time -v  $MINIMAP -a -x $PROFILE hg38noAlt_$i.idx $READS  --multi-prefix tmp > merge.$i.sam 2>> $SAM_MAP_LOG
		echo "Mapping to partitioned index  for $i parts with merging: PAF"
		/usr/bin/time -v $MINIMAP -x $PROFILE hg38noAlt_$i.idx $READS  --multi-prefix tmp > merge.$i.paf 2>> $PAF_MAP_LOG
	done

	rm *.sam *.paf *.idx
}


logcat(){

	echo -e "type\t1 part\t2 part\t4 part\t8 part\t 16 part"
	echo -e -n "\nindex_ont\t"
	cat index_ont.log | grep "Maximum resident" | awk '{print $NF/(1024*1024)}' | tr '\n' '\t'
	echo -e -n "\nsam_ont\t"
	cat sam_ont.log | grep "Maximum resident" | awk '{print $NF/(1024*1024)}' | tr '\n' '\t'
	echo -e -n "\npaf_ont\t"
	cat paf_ont.log | grep "Maximum resident" | awk '{print $NF/(1024*1024)}' | tr '\n' '\t'
	echo ""
}


stats_after_opti "ont"
logcat