#!/bin/bash

#compile for the exact match based compare mode
g++ -Wall -O2 -std=c++11  comparesam.c -o exact_compare_primary_only.out

#compile for the overlap based compare mode
g++ -Wall -O2 -std=c++11 -DOVERLAP_BASED_EVAL=1 comparesam.c -o overlap_compare_primary_only.out

UNIPART=unipart.sam
PARTIDX=multipart_merge.sam

mkdir exact_compare overlap_compare

./exact_compare_primary_only.out $UNIPART $PARTIDX > exact_compare/summary.txt
mv mismatches.tsv mismatches_a.bed mismatches_b.bed only_in_a.bed only_in_a.tsv only_in_b.bed only_in_b.tsv exact_compare

./overlap_compare_primary_only.out $UNIPART $PARTIDX > overlap_compare/summary.txt
mv mismatches.tsv mismatches_a.bed mismatches_b.bed only_in_a.bed only_in_a.tsv only_in_b.bed only_in_b.tsv overlap_compare
