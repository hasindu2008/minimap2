#!/bin/bash

#by default 0. If you want to get peak ram usage make this 1. In that case make sure /usr/bin/time is installed.
MEASURE_MAX_RAM=0;

error_exit ()
{
	echo $1
	exit 1
}


if [ "$#" -ne 5 ]; then
    echo "usage : $0 <reference.fa> <num_parts> <out.idx> <minimap2_exe> <minimap2_profile>"
	echo ""
	echo "reference.fa - path to the fasta file containing the reference genome"
	echo "num_parts - number of partitions in the index"
	echo "out.idx - path to the file to which the index should be dumped"
	echo "minimap2_exe - path to the minimap2 executable"
	echo "minimap2_profile - minimap2 pre-set for indexing (map-pb or map-ont)"
	echo ""
	echo "Example : $0 hg19.fa 4 hg19.idx minimap2 map-ont"
    exit 1
fi

REFERENCE=$1
NUMPARTS=$2
OUTIDX=$3
MINIMAP=$4
PROFILE=$5

#checking args
[ -e "$REFERENCE" ] || error_exit "ERROR : $REFERENCE does not exist."
[ "$NUMPARTS" -lt 2 ] && error_exit "num_parts should be 2 or higher"
command "$MINIMAP" --version >/dev/null 2>&1 || error_exit "ERROR : $MINIMAP does not exist or does not have executable permission"
[ "$PROFILE" = "map-ont" ] || [ "$PROFILE" = "map-pb" ]  || error_exit "ERROR : minimap2_profile accepts only map-ont or map-pb"
[ -w "./" ] || error_exit "ERROR : Current directory is not writable."


BASEDIR=$(dirname "$0")

set -e

echo "Compiling divide.c"
gcc -Wall -O2 $BASEDIR/divide.c -o divide

echo "Running divider"
if [ "$MEASURE_MAX_RAM" -eq 1 ]; then
	/usr/bin/time -v ./divide $REFERENCE $NUMPARTS  2> >(tee divider.log)
else
	./divide $REFERENCE $NUMPARTS
fi

echo "Minimap2 Indexing"
LIST="";
NUMPARTS2=$(echo "$NUMPARTS-1" | bc )
for i in $(seq 0 $NUMPARTS2)
do
	echo "Creating partition $i"
	if [ "$MEASURE_MAX_RAM" -eq 1 ]; then
		/usr/bin/time -v $MINIMAP -a -x $PROFILE -d part$i.idx part$i.fa 2> >(tee part$i.log)
	else
		$MINIMAP -a -x $PROFILE -d part$i.idx part$i.fa > /dev/null
	fi
	LIST="$LIST"" part$i.idx"
	rm part$i.fa
done


echo "concatenating $LIST"
cat $LIST > $OUTIDX
rm $LIST

if [ "$MEASURE_MAX_RAM" -eq 1 ]; then
	max_mem=$(grep "Maximum resident set size" *.log | awk '{print $NF/1024}' | sort -n | tail -1);
	echo "Peak memory usage : $max_mem MB"
	for i in $(seq 0 $NUMPARTS2)
	do
		rm part$i.log
	done
fi
