# Evaluation scripts for partitioned indexes

## `synth_eval.sh`

The bash shell script for evaluating the mapping accuracy from the synthetic reads.
Performs the mapping using minimap2 and evaluation using Paftools for
single-idx, part-idx-merged and part-idx-no-merge.
See the comments in  the script for more information.

### `synth_plot.m`
The Matlab script that generate the plots the fraction of mapped reads against the error rate of mapped reads.
Should be run after `synth_eval.sh`. Tested on Matlab 2017.

## `chimeric_eval.sh`
The shell script for mapping the chromothriptic read with different options.
Performs the mapping using minimap2 for single-idx, part-idx-merged and part-idx-no-merge.
See the comments in  the script for more information.

### `chimeric_plot.m`
The Matlab script that generates the plots to visualise the mappings of a chromothriptic read against the reference.
Should be run after `chimeric_eval.sh`. Tested on Matlab 2017.

## `real_get_stat.sh`
The shell script for mapping a set of real reads and deriving mapping statistics.
Performs the mapping using minimap2 for single-idx, part-idx-merged and part-idx-no-merge.
See the comments in  the script for more information.

## `run_comparesam.sh`
The shell script to compile and run `comparesam.c` to thoroughly evaluate the difference between the mapping outputs from single-idx and part-idx-merged. Useful when no truth set is available (for instance for real datasets), but only required to compare the difference between the output from two mappers.

### `comparesam.c`
The C program that compares two sam files.
Let  *a* and *b* be two samfiles for the same set of reads,
but mapped with different mappers (or same mapper with different options).
comparesam will compare *a* and *b* and give mappings statistics such as the number of reads which are :
- unmapped in both
- correct - the read maps to the same location (locations overlaps if overlap based evaluation is set)
- incorrect - the read DO NOT map to the same location (locations DO NOT overlap if overlap based evaluation is set)
- only mapped in  *a*
- only mapped in  *b*
- primary mapping in *a* is a supplementary mapping in *b*
- primary mapping in *b* is a secondary mappings in *b*

comparesam will also output (as tsv files and bed file) reads that :
- mismatch between *a* and *b*
- unique to *a*
- unique to *b*

functionality is as follows :
- samfile *a* and *b* are sequentially read while loading all the mappings for a
particular read name at a time.
- IMPORTANT : the samfiles *a* and *b* should have the reads in the same order and the multiple mappings for a given read ID should be adjacently located.
- For each loaded read, we compare the mappings between *a* and *b*. In *a* always the primary mapping is considered despite the value of
CONSIDER_SUPPLEMENTARY and CONSIDER_SECONDARY flags.
- In *b* supplementary and secondary mappings can also be considered if the flags CONSIDER_SUPPLEMENTARY and CONSIDER_SECONDARY are set.
- The comparison statistics are updated during the comparison and
the entries will be written to the tsv and bed file if required.
- Finally we print the comparison statistics.

## `mem_usage_no_chr_balance.sh`
The bash script that evaluates the peak memory usage for an index without chromosome balancing.
Uses the inbuilt indexing option in minimap2

## `Mem_usage_with_chr_balance.sh`
The bash script that evaluates the peak memory usage for the index with chromosome balancing.
Uses `divide_and_index.sh` to construct the index.

## `gettime.sh`
The bash script that evaluates the runtime of single-idx and part-idx-merged.
