#!/bin/bash

JEME_DIR=$1
INPUT_BED=$2

mkdir overlaps/JEME

for F in $JEME_DIR/*.csv.gz; do
  fname=`basename $F`
  ID=$(echo $fname | cut -f 2 -d '.')
  bedtools intersect -wb -a $INPUT_BED -b <(zcat $F | perl -F',' -ane '@x=split(/:|-/, $F[0]); print join("\t", @x, $F[1], $F[2]);') \
    | cut -f 1,3,8,9 > overlaps/JEME/$ID.overlap.txt
done

