#!/bin/bash


#path to the modified minimap2 binary
MINIMAP=../../minimap2-arm/minimap2

#path to chromosome balancing indexer c file
DIVIDER=../../minimap2-arm/misc/idxtools/divide.c


#path to the reference index (download and extract from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz )
REFERENCE=hg38noAlt.fa


#read file (download and extract from http://s3.amazonaws.com/nanopore-human-wgs/rel3-nanopore-wgs-84868110-FAF01132.fastq.gz )
READS=rel3-nanopore-wgs-84868110-FAF01132.fastq



#############################################################################

PROFILE=map-ont


set -e


#generate indexes
GENERATE_IDX(){


	echo "Compiling divide.c"
	gcc -Wall -O2 $DIVIDER -o divide 

	#INDEXING for single-idx
	echo "Generating index"
	/usr/bin/time -v $MINIMAP -a -x $PROFILE $REFERENCE -d hg38noAlt_single.idx > /dev/null 2> index_1.log

	#indexing for n-part-idx	
	for i in 2 4 8 16
	do
		echo "Generating index for $i"
		echo "Running divider"
		/usr/bin/time -v ./divide $REFERENCE $i  2> index_$i.log

		echo "Minimap2 Indexing"
		LIST="";
		NUMPARTS2=$(echo "$i-1" | bc )
		for j in $(seq 0 $NUMPARTS2)
		do
			echo "Creating partition $j"
			/usr/bin/time -v $MINIMAP -a -x $PROFILE -d part$j.idx part$j.fa > /dev/null 2>> index_$i.log
			LIST="$LIST"" part$j.idx"
			rm part$j.fa
		done

		echo "concatenating $LIST"	
		/usr/bin/time -v cat $LIST > hg38noAlt_$i.idx 2>> index_$i.log
		rm $LIST

	done

}

	

#SAM MAPPING
MAP(){

	SAM_MAP_LOG=map.log	
	echo "Sam map log" > $SAM_MAP_LOG
	
	#MAPPING with single-idx
	echo "Generating SAM"
	/usr/bin/time -v $MINIMAP -a -x $PROFILE -t8 hg38noAlt_single.idx $READS  > groundtruth.sam 2>> $SAM_MAP_LOG
	
	#MAPPING for n-part-idx
	for i in 2 4 8 16
	do
		echo "Multi-part index with modified minimap2 for $i with merging: SAM"
		/usr/bin/time -v  $MINIMAP -a -x $PROFILE -t8 hg38noAlt_$i.idx $READS  --multi-prefix tmp > merge.$i.sam 2>> $SAM_MAP_LOG
	done

	#rm *.idx *.csv *.sam
	rm *.csv *.sam
}


#PAF MAPPING
MAP_PAF(){

	PAF_MAP_LOG=map_paf.log	
	echo "paf map log" > $PAF_MAP_LOG

	#MAPPING
	echo "Generating PAF"
	/usr/bin/time -v $MINIMAP -x $PROFILE -t8 hg38noAlt_single.idx $READS  > groundtruth.paf 2>> $PAF_MAP_LOG
	
	for i in 2 4 8 16
	do
		echo "Multi-part index with modified minimap2 for $i with merging: PAF"
		/usr/bin/time -v  $MINIMAP -x $PROFILE -t8 hg38noAlt_$i.idx $READS  --multi-prefix tmp > merge.$i.paf 2>> $PAF_MAP_LOG
	done

	rm *.idx  *.paf

}

#get iformation from logs
logcat(){
	
	echo -e -n "indexing_1_part\t"
	grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" index_1.log | cut -d ' ' -f 8 |tr ':' \\t |  awk '{if(NF==1) print; else{ if(NF==2) print(($1*60)+$2); else print(($1*3600)+($2*60)+$3)}}' | tr '\n' '\t'

	for i in 2 4 8 16
	do
		echo -e -n "\nindexing_$i""_parts\t"
		grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" index_$i.log | cut -d ' ' -f 8 |tr ':' \\t |  awk '{if(NF==1) print; else{ if(NF==2) print(($1*60)+$2); else print(($1*3600)+($2*60)+$3)}}' | tr '\n' '\t'
	done

	echo -e -n "\n\ntime for mapping SAM (merge time not included)(1,2,4,8,16parts)\t"
	grep -B3 "Real time" map.log | grep "worker_pipeline" | awk -F'[:*]' '{print $5}' |  tr '\n' '\t'

	echo -e -n "\ntime for mapping with merge time SAM (1,2,4,8,16parts)\t"
	grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" map.log | cut -d ' ' -f 8 |tr ':' \\t |  awk '{if(NF==1) print; else{ if(NF==2) print(($1*60)+$2); else print(($1*3600)+($2*60)+$3)}}' | tr '\n' '\t'

	echo -e -n "\n\ntime for mapping PAF (merge time not included)(1,2,4,8,16parts)\t"
	grep -B3 "Real time" map_paf.log | grep "worker_pipeline" | awk -F'[:*]' '{print $5}' |  tr '\n' '\t'

	echo -e -n "\ntime for mapping with merge time PAF (1,2,4,8,16parts)\t"
	grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" map_paf.log | cut -d ' ' -f 8 |tr ':' \\t |  awk '{if(NF==1) print; else{ if(NF==2) print(($1*60)+$2); else print(($1*3600)+($2*60)+$3)}}' | tr '\n' '\t'


	echo ""
}


#get additional information from logs
logcat2(){

	echo "sam idxload_abs"
	grep "M::mm_idx_stat::" map.log  |   awk -F'[:*]' 'BEGIN{i=1}{if(i==1||i==1+2||i==1+2+4||i==1+2+4+8||i==1+2+4+8+16){printf $5"\n"}else{printf $5"\t"};i=i+1}'
	echo "sam worker_pipeline_prev_abs"	
	grep "M::mm_idx_stat::" map.log -B4 | grep "worker_pipeline" |   awk -F'[:*]' 'BEGIN{i=1;printf "0\n0\t"}{if(i==1||i==1+3||i==1+3+7){printf $5"\n0\t"}else{printf $5"\t"};i=i+1}'
	echo ""	
	echo "paf idxload_abs"
	grep "M::mm_idx_stat::" map_paf.log  |   awk -F'[:*]' 'BEGIN{i=1}{if(i==1||i==1+2||i==1+2+4||i==1+2+4+8||i==1+2+4+8+16){printf $5"\n"}else{printf $5"\t"};i=i+1}'
	echo "paf worker_pipeline_prev_abs"
	grep "M::mm_idx_stat::" map_paf.log -B4 | grep "worker_pipeline" |   awk -F'[:*]' 'BEGIN{i=1;printf "0\n0\t"}{if(i==1||i==1+3||i==1+3+7){printf $5"\n0\t"}else{printf $5"\t"};i=i+1}'
	
}


GENERATE_IDX
echo "Indexes generated"
MAP
echo "Mapping sam done"
MAP_PAF
echo "Mapping paf done"
logcat
logcat2



