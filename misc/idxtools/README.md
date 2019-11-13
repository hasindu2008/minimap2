# Index construction with chromosome size balancing

`divide_and_index.sh` is the wrapper script for balanced index construction.
It takes the reference genome and outputs a partitioned index optimised for reduced peak memory.

Its usage is as follows.

```
usage : ./divide_and_index.sh <reference.fa> <num_parts> <out.idx> <minimap2_exe> <minimap2_profile>

reference.fa - path to the fasta file containing the reference genome
num_parts - number of partitions in the index
out.idx - path to the file to which the index should be dumped
minimap2_exe - path to the minimap2 executable
minimap2_profile - minimap2 pre-set for indexing (map-pb or map-ont)

Example : ./divide_and_index.sh hg19.fa 4 hg19.idx minimap2 map-ont

```

Functionality of `divide_and_index.sh` is as follows.
1. Compiling `divide.c` using `gcc` to produce `divide`.
2. Call the compiled binary `divide` to split the reference genome into partitions such that the total length of chromosomes in each partition are approximately equal.
3. Calling the minimap2 binary separately on each reference partition to produce a separate index file for each partition.
4. Combining all the index files to produce a single partitioned index file.


# Running Minimap2 on a partitioned index with merging

To run minimap2 on an index created using the above method :
```
minimap2 -x <profile> <partioned_index.idx> <reads.fastq> --multi-prefix <tmp-prefix>
```
`--multi-prefix` which takes a prefix for temporary files, enables the merging of the outputs generated through iterative mapping to index partitions.


# Example

1. Download and compile minimap2 that supports partitioned indexes and merging
```
git clone https://github.com/hasindu2008/minimap2-arm && cd minimap2-arm && make
```

1. Download the human reference genome and create a partitioned index with 4 partitions
```
wget -O hg38noAlt.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz && gunzip hg38noAlt.fa.gz
./misc/idxtools/divide_and_index.sh hg38noAlt.fa 4 hg38noAlt.idx ./minimap2 map-ont
```

2. Download a Nanopore NA12878 dataset and run Minimap2 with merging
```
wget http://s3.amazonaws.com/nanopore-human-wgs/rel3-nanopore-wgs-84868110-FAF01132.fastq.gz
./minimap2 -a -x map-ont hg38noAlt.idx rel3-nanopore-wgs-84868110-FAF01132.fastq.gz --multi-prefix tmp > aligned_out.sam
```

Notes :

- To perform mapping without base-level alignment use:
```
./minimap2 -x map-ont hg38noAlt.idx rel3-nanopore-wgs-84868110-FAF01132.fastq.gz --multi-prefix tmp > aligned_out.paf
```
- From [Minimap2 version 2.12-r827](https://github.com/lh3/minimap2/blob/master/NEWS.md#release-212-r827-6-august-2018) onwards, the merging functionality has been integrated into the main repository. This version additionally supports paired-end short reads and the merging operation is multi-threaded. Use `--split-prefix` option instead of `--multi-prefix`.

- Some useful evaluation scripts (accuracy, difference, memory usage and runtime) are [here](eval)
