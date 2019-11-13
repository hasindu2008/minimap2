#!/bin/bash


#path to the modified minimap2 binary
MINIMAP=../minimap2-arm/minimap2

#path to the reference index (download and extract from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz )
REFERENCE=hg38noAlt.fa

#read file (download from http://s3.amazonaws.com/nanopore-human-wgs/rel3-nanopore-wgs-84868110-FAF01132.fastq.gz )
fastq=rel3-nanopore-wgs-84868110-FAF01132.fastq.gz



#############################################################################


file1=${fastq%%.*}
name=$(basename $file1)

i=140M #16 parts in the index

SIMPLE_STATS(){
	echo "Num bases"
	zcat $fastq | awk 'BEGIN{sum=0}{if(NR%4==2) {sum=sum+length($0); }}END{print sum}'
	echo ""
	echo "Num reads"
	echo "$(cat $fastq | wc -l)/4" | bc 
	echo ""
}


GENERATE_IDX(){

#generate indices
echo "Generating single reference index"
$MINIMAP -a -x map-ont $REFERENCE -d hg38noAlt_single.idx > /dev/null 2> hg38noAlt_single.idx.log

echo "Generating partitioned index with 16 parts"
$MINIMAP -a -x map-ont $REFERENCE -d hg38noAlt_$i.idx -I$i > /dev/null 2> hg38noAlt_$i.idx.log

}

	
RUN_UNIPART(){
	/usr/bin/time -v $MINIMAP -a -x map-ont hg38noAlt_single.idx $fastq -t 8 > unipart.sam 2> unipart.sam.log
}

MULTIPART_NOMERGE(){
	/usr/bin/time -v $MINIMAP -a -x map-ont hg38noAlt_$i.idx $fastq -t 8 > before.sam 2> before.sam.log
}


MULTIPART_MERGE(){
	/usr/bin/time -v $MINIMAP -a -x map-ont hg38noAlt_$i.idx $fastq -t 8 --multi-prefix tmp > multipart_merge.sam 2> multipart_merge.sam.log
}

MORE_STATS(){
	for file in unipart.sam multipart_merge.sam
	do
		echo "File $file"
		echo "file size"
		ls -lh $file | awk '{print $5}'
		echo "num entries"
		samtools view $file | wc -l
		echo "flag stats"
		samtools view $file | awk '{print $2}' | sort -n | uniq -c | awk '{print $2,$1}'
		echo "qual_scores : all except unmapped"
		samtools view -F4 $file |  awk '{print $5}' | sort -n | uniq -c | awk '{print $2,$1}'
		echo "qual scores : primary"
		samtools view -F260 $file |  awk '{print $5}' | sort -n | uniq -c | awk '{print $2,$1}'
		echo ""
		echo ""
	done

	file=before.sam
	echo "File $file"
	echo "file size"
	ls -lh $file | awk '{print $5}'
	echo "num entries"
	cat hg38noAlt.header.sam $file | samtools view - | wc -l
	echo "flag stats"
	cat hg38noAlt.header.sam $file | samtools view - | awk '{print $2}' | sort -n | uniq -c | awk '{print $2,$1}'
	echo "qual_scores : all except unmapped"
	cat hg38noAlt.header.sam $file |  samtools view -F4 - |  awk '{print $5}' | sort -n | uniq -c | awk '{print $2,$1}'
	echo "qual scores : primary"
	cat hg38noAlt.header.sam $file | samtools view -F260 - |  awk '{print $5}' | sort -n | uniq -c | awk '{print $2,$1}'
	echo ""
	echo ""
}



SIMPLE_STATS
echo "Simple stats generated"
GENERATE_IDX
echo "Indexes generated"
RUN_UNIPART
echo "Single reference index done"
MULTIPART_NOMERGE
echo "Partitioned index without merge done"
MULTIPART_MERGE
echo "Partitioned index with merge done"
MORE_STATS > sam_stat.txt
echo "Stats generated"

